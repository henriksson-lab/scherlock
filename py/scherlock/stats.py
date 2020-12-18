from pandasql import sqldf

def cellcounts(adata, columns):
    """Get the number of cells for each of the condition.

    Args:
        adata (anndata): The single-cell count object
        columns (list of strings): The columns specifying conditions

    Returns:
        Dataframe with counts for each condition. The conditions are sorted according to
        their order in the columns argument.
    """
    tobs=adata.obs
    jcol = ",".join(columns)
    return(sqldf("select " + jcol+ ",count(*) as count from tobs group by "+jcol+" order by "+jcol, locals()))

#cellcounts(adata, ["genotype", "treatment"])

