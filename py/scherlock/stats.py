from pandasql import sqldf

def cellcounts(adata, columns, groupby=[], percentage=False):
    """Get the number of cells for each of the condition.

    Args:
        adata (anndata): The single-cell count object
        columns (list of strings): The columns specifying conditions
        percentage (boolean): Calculate percentage rather than absolute counts
        groupby (list of strings): If percentage, calculate percentage of cells within these groups

    Returns:
        Dataframe with counts for each condition. The conditions are sorted according to
        their order in the columns argument.
    """
    tobs=adata.obs
    jcol = ",".join(columns)
    if not percentage:
        ### Just show a raw count
        return(sqldf("select " + jcol+ ",count(*) as count from tobs group by "+jcol+" order by "+jcol, locals()))
    else:
        ### Show percentage
        jcol_count = ",".join(columns+["count"])

        #Total number of cells in the group
        if len(groupby)==0:
            df_group = sqldf("select count(*) as totcount from tobs", locals())
        else:
            gcol = ",".join(groupby)
            df_group = sqldf("select " + gcol + ",count(*) as totcount from tobs group by "+gcol, locals())

        #Total number of cells in each condition
        df_cond = sqldf("select " + jcol+ ",count(*) as condcount from tobs group by "+jcol+" order by "+jcol, locals())

        #Compare cell counts
        if len(groupby)==0:
            df = df_cond
            df["totcount"] = tobs.shape[0]
        else:
            df = df_group.merge(df_cond)
        jcol_perc = ",".join(columns+["cast (condcount as real)/totcount as perc"])
        df = sqldf("select " + jcol_perc+ " from df order by "+jcol, locals())
        return(df)


### Calculate absolute cell counts
#cellcounts(adata, ["genotype", "treatment"])

### Calculate percentage cells of the whole
#cellcounts(adata, ["genotype", "treatment"], percentage=True)

### Calculate percentage cells, within each group
cellcounts(adata, ["genotype", "treatment"], percentage=True, groupby=["genotype"])


