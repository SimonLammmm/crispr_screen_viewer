from pathlib import Path
import re
import pandas as pd
import os, pickle
import sqlite3
from typing import Dict, Collection
from crispr_screen_viewer.functions_etc import (
    timepoint_labels,
    doi_to_link,
    index_of_true,
    LOG,
)

class DataSet:
    """Class for holding, retrieving screen data and metadata.

    Storage and retreval of amalgamated results tables. Experiment tables as
    single analysis type/statistic with filename {ans_type}_{stat}.csv
    Currently supports MAGeCK ('mag') and DrugZ ('drz').

    Attributes:
        exp_data: data tables keyed first by analysis type then 'score'|'fdr'
        comparisons: descriptions of exp_data samples. Indexed by comparison keys
        score/analysis_labels: text labels for supported analysis types
        genes: Index used in tables in exp_data

    Methods:
        get_results: get dict of score and fdr for specified analysis types and
            datasets."""
    def __init__(self, source_directory, print_validations=True):
        LOG.debug(source_directory)
        source_directory = Path(source_directory)
        # for future ref
        self.source_directory = str(source_directory)

        # shorthand internal name: label name
        avail_analyses = []
        for ans in ('drz', 'mag'):
             #if os.path.isfile(source_directory/f"{ans}_fdr.csv.gz"): # no longer using
             avail_analyses.append(ans)

        self.available_analyses = avail_analyses
        self.analysis_labels = {'drz':'DrugZ', 'mag':'MAGeCK'}
        self.score_labels = {'mag':'Log2(FC)', 'drz':'NormZ'}

        # # put the data tables in {analysis_type:{score/fdr:pd.DataFrame}} format dictionary
        # exp_data = {ans:{stt:pd.read_csv(source_directory/f"{ans}_{stt}.csv.gz", index_col=0)
        #                  for stt in ('score', 'fdr')}
        #             for ans in self.available_analyses}

        # unify the indexes
        
        
        genes = pd.Index([])
        # use an index from a table from each analysis
        for analysis in self.available_analyses:
        # genes = genes.union(exp_data[analysis]['fdr'].index)
            self.con = sqlite3.connect(self.source_directory + '/ddrcs.db', uri=True)
            self.c = self.con.cursor()
            genes_cursor = self.c.execute("SELECT [gene] FROM " + analysis + "_score")
            genes1 = genes_cursor.fetchall()
            self.con.close()
            genes0 = [str(g[0]) for g in genes1]
            genes = genes.union(genes0)
        # remove non-targeting and sequences from the genes list
        r = re.compile("^(?!.*Non-targeting.*)")
        genes = list(filter(r.match, list(genes)))
        r = re.compile("^(?!.*[ATCG]{19}.*)")
        genes = list(filter(r.match, list(genes)))
        self.genes = pd.Index(genes)
        # # reindex with the union of genes
        # self.exp_data = {ans:{stt:exp_data[ans][stt].reindex(genes)
        #                     for stt in ('score', 'fdr')}
        #                for ans in self.available_analyses}

        comparisons = pd.read_csv(source_directory/'comparisons_metadata.csv.gz', )
        # this is sometimes put in wrong...
        m = comparisons['Timepoint'] == 'endpoint'
        comparisons.loc[m, 'Timepoint'] = 'endpoints'
        comparisons.loc[m, 'Control group'] = comparisons.loc[m, 'Control group'] \
            .apply(lambda x: x.replace('endpoint', 'endpoints'))

        comparisons = comparisons.set_index('Comparison ID', drop=False)
        #comparisons.loc[:, 'Available analyses'] = comparisons['Available analyses'].str.split('|')
        try:
            comparisons = comparisons.drop('Available analyses', axis=1)
        except KeyError:
            pass

        # fill in blank treatments
        comparisons.loc[comparisons.Treatment.isna(), 'Treatment'] = 'No treatment'
        # these cols could be blank and aren't essential to have values
        for col in  ['Cell line', 'Library', 'Source']:
            if col not in comparisons.columns:
                comparisons.loc[:, col] = 'Unspecified'
            comparisons.loc[comparisons[col].isna(), col] = 'Unspecified'

        # replace some values with ones that read better
        for old, new in timepoint_labels.items():
            comparisons.loc[comparisons.Timepoint == old, 'Timepoint'] = new

        # list of all datasources for filtering
        self.data_sources = comparisons.Source.fillna('Unspecified').unique()
        # main metadata tables
        self.comparisons = comparisons
        self.experiments_metadata = pd.read_csv(f'{source_directory}/experiments_metadata.csv.gz', )
        # rename "Experiment name" to "Experiment ID" for consistency
        colmap = {k:k for k in self.experiments_metadata}
        colmap['Experiment name'] = 'Experiment ID'
        self.experiments_metadata.columns = self.experiments_metadata.columns.map(colmap)
        self.experiments_metadata.set_index('Experiment ID', drop=False, inplace=True)

        # add formated DOI to the comparisons metadata
        dois = self.experiments_metadata.loc[
            self.comparisons['Experiment ID'],
            'DOI'
        ].apply(doi_to_link).values
        self.comparisons.insert(2, 'DOI', dois)

        # add citation to comparisons table
        try:
            cites = self.experiments_metadata.loc[
                self.comparisons['Experiment ID'],
                'Citation'
            ].values
            self.comparisons.loc[:, 'Citation'] =  cites
        except:
            LOG.warning('Citations column missing from exeriments_metadata')
            self.comparisons.loc[:, 'Citation'] = ''


        # DF of previous symbols and IDs for currently used.
        try:
            pidf = pd.read_csv(
                os.path.join(source_directory, 'previous_and_id.csv.gz'), index_col=0
            )
            # switching to using the default names in hgnc_complete_text.txt
            #   so here's another inconsistancy fix...
            # (let's hope HUGO never changes the table format hah)
            cmap = {k:k for k in pidf.columns}
            cmap['HGNC_ID'] = 'hgnc_id'
            cmap['Previous_symbol'] = 'prev_symbol'
            pidf.columns = pidf.columns.map(cmap)

            self.previous_and_id = pidf.fillna('')

        except FileNotFoundError:
            LOG.warning("file 'previous_and_id.csv.gz' is missing.")
            # when .loc fails to find a name in the table it just uses the current name.
            self.previous_and_id = pd.DataFrame()

        self.previous_and_id.fillna('', inplace=True)

        if print_validations:
            self.validate_comparisons()
            self.validate_previous_and_id()

    def validate_comparisons(self):
        """Print information that might be helpful in spotting data validity issues
        Check for comparisons present in the metadata/actual-data but missing
          in the other"""
        all_good = True
        for ans in self.available_analyses:
            self.con = sqlite3.connect(self.source_directory + '/ddrcs.db', uri=True)
            self.c = self.con.cursor()
            score_comps_cursor = self.c.execute("SELECT * FROM " + ans + "_score")
            score_comps = [d[0] for d in score_comps_cursor.description]
            self.con.close()
            score_comps = pd.DataFrame(index = score_comps).index
            #score_comps = self.exp_data[ans]['score'].columns
            meta_comps = self.comparisons.index

            meta_in_score = meta_comps.isin(score_comps)
            missing_in_data = meta_comps[~meta_in_score]
            # todo log.warning
            # todo check experiments metadata
            if missing_in_data.shape[0] > 0:
                all_good = False
                print(
                    f"Comparisons in comparisons metadata but not in {ans}_score.csv.gz:"
                    f"\n    {', '.join(missing_in_data)}\n"
                )
            score_in_meta = score_comps.isin(meta_comps)
            missing_in_score = score_comps[~score_in_meta]
            if missing_in_data.shape[0] > 0:
                all_good = False
                print(
                    f"Comparisons in {ans}_score.csv.gz, but not in comparisons metadata:"
                    f"\n    {', '.join(missing_in_score)}\n"
                )
        comps = self.comparisons.index
        if comps.duplicated().any():
            all_good = False
            print('Duplicate comparisons found - this will probably stop the server from working:')
            print('   ' , ', '.join(sorted(comps.index[comps.index.duplicated(keep=False)])))
        if all_good:
            print(f'All comparisons data in {self.source_directory} are consistent')

        # check all comparison ExpID appear in experiments metadata
        # the reverse isn't fatal
        expids = self.comparisons['Experiment ID'].unique()
        found = [xi in self.experiments_metadata.index for xi in expids]
        if not all(found):
            not_found = [x for (x, b) in zip(expids, found) if not b]
            print('Experiment IDs used in comparisons_metadata not found in experiments_metadata:\n'
                  f'   {", ".join(not_found)}')

    def validate_previous_and_id(self):
        """Check all stats index in datasets are in previous_and_id"""
        for ans in self.available_analyses:
            self.con = sqlite3.connect(self.source_directory + '/ddrcs.db', uri=True)
            self.c = self.con.cursor()
            score_genes_cursor = self.c.execute("SELECT [gene] FROM " + ans + "_score")
            score_genes = score_genes_cursor.fetchall()
            self.con.close()
            score_genes0 = [g[0] for g in score_genes]
            score_index = pd.DataFrame(score_genes0, index = score_genes0)
            #score_index = self.exp_data[ans]['score'].index
            m = score_index.index.isin(self.previous_and_id.index)
            print(m.sum(), 'of', len(m), f'gene symbols have record in previous_and_id.csv.gz, in file {ans}_score.csv.gz')


    def get_score_fdr(self, score_anls:str, fdr_anls:str=None,
                      comparisons:list=None, genes:list=None,
                      data_sources:Collection= 'all') -> Dict[str, pd.DataFrame]:
        """Get score and FDR tables for the analysis types & data sets
        for the genes and comparisons specified.
        Tables give the per gene values for included comparisons.

        Arguments:
            score_anls: The analysis type from which to get the score values
                per gene
            fdr_anls: Optional. As score_anslys
            comparisons: Optional. Return only the comparisons specified.
                Default None. If None, then all comparisons are returned.
            genes: Optional. Return only the genes specified. Default None.
                If None, then return all genes.
            data_sources: Data sources (i.e. SPJ, or other peoples papers) to
                include in the returned DFs. Any comparison that comes from a
                dataset that does not have both fdr/score analysis types
                available will not be present in the table.

        Returns {'score':pd.DataFrame, 'fdr':pd.DataFrame}"""


        # todo surely we don't need to do this every time?
        #   write tables for each score/stat, store in a dict
        
        # if only one type supplied, copy it across
        if fdr_anls is None:
            fdr_anls = score_anls

        # SQL query builder: score
        query_score = "SELECT "
        # comparisons control
        if comparisons is None:
            query_score = query_score + "*"
        else:
            query_score_comparisons = comparisons.copy()
            query_score_comparisons.insert(0, "gene") 
            query_score = query_score + "[" + "],[".join(query_score_comparisons) + "]" 
        # method control
        query_score_anls = score_anls + "_score"
        query_score = query_score + " FROM " + query_score_anls
        # genes control
        self.con = sqlite3.connect(self.source_directory + '/ddrcs.db', uri=True)
        self.c = self.con.cursor()
        if genes is not None:
            query_score = query_score + " WHERE gene = ?"
            for x in range(1, len(genes)):
                query_score = query_score + " OR gene = ?"
            scoreCursor = self.c.execute(query_score, (genes))
        else:
            scoreCursor = self.c.execute(query_score)
                
        score = pd.DataFrame(scoreCursor.fetchall(), columns = [d[0] for d in scoreCursor.description])
        self.con.close()
        score = score.set_index('gene')
        
        # SQL query builder: fdr
        query_fdr = "SELECT "
        # comparisons control
        if comparisons is None:
            query_fdr = query_fdr + "*"
        else:
            query_fdr_comparisons = comparisons.copy()
            query_fdr_comparisons.insert(0, "gene") 
            query_fdr = query_fdr + "[" + "],[".join(query_fdr_comparisons) + "]" 
        # method control
        query_fdr_anls = fdr_anls + "_fdr"
        query_fdr = query_fdr + " FROM " + query_fdr_anls
        # genes control
        self.con = sqlite3.connect(self.source_directory + '/ddrcs.db', uri=True)
        self.c = self.con.cursor()
        if genes is not None:
            query_fdr = query_fdr + " WHERE gene = ?"
            for x in range(1, len(genes)):
                query_fdr = query_fdr + " OR gene = ?"
            fdrCursor = self.c.execute(query_fdr, (genes))
        else:
            fdrCursor = self.c.execute(query_fdr)
            
        fdr = pd.DataFrame(fdrCursor.fetchall(), columns = [d[0] for d in fdrCursor.description])
        self.con.close()
        fdr = fdr.set_index('gene')
        
        score_fdr = {'score':score, 'fdr':fdr}
        
        # score_fdr = {stt:self.exp_data[ans][stt] for ans, stt in ((score_anls, 'score'), (fdr_anls, 'fdr'))}     

        if data_sources == 'all':
            return score_fdr

        # Filter returned comparisons (columns) by inclusion in data sources and having
        #   results for both analysis types
        comps_mask = self.comparisons.Source.isin(data_sources)
        # for analysis_type in (score_anls, fdr_anls):
        #     m = self.comparisons['Available analyses'].apply(lambda available: analysis_type in available)
        #     comps_mask = comps_mask & m
        comparisons = index_of_true(comps_mask)
        score_fdr = {k:tab.reindex(columns=comparisons) for k, tab in score_fdr.items()}

        return score_fdr

    def dropdown_gene_label(self, gn):
        try:
            row = self.previous_and_id.loc[gn]
        except:
            return gn

        if not row.hgnc_id:
            return gn

        s = f"{gn}  ({row.hgnc_id}"
        s_end = ')'
        if row.prev_symbol:
            s_end = f"; {row.prev_symbol})"

        return s + s_end

def load_dataset(paff):
    """If paff is a dir, the dataset is constructed from the files
    within, otherwise it is assumed to be a pickle."""
    if os.path.isfile(paff):
        LOG.info('args.data_path is a file, assuming pickle and loading.')
        with open(paff, 'rb') as f:
            data_set = pickle.load(f)
    else:
        data_set = DataSet(Path(paff))

    return data_set
