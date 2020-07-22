def _dropped(reason, clusterIDs, log_path="logs/dropped.txt"):
    """ tracked dropped clusterIDs and rationale to a single log file 
    Manually curated by pipeline runner
    Downstream steps will drop these clusterIDs from analysis
    
    """
    with open(log_path, "a+") as f:
        for _id in clusterIDs:
            f.write(f"{_id},{reason}")
