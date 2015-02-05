Right now, the script for running a simplified version of this is "Fiber_Analyzer_Run.m", while the older version is "FADTTS_analyze.m" However, these are designed to be able to run independently. Sample group and profile files are also provided.

The following calls can be made to individual functions in this library:

    [something]_plot('/path/to/FA_file.csv','/path/to/groups_file.csv', 'Fiber Name')

This is the easiest way to compute local and global p-values, averages and standard deviations, and raw values, for any given set of fiber data.
