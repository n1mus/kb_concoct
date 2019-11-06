/*
A KBase module: kb_concoct
*/

module kb_concoct {

    /* A boolean - 0 for false, 1 for true.
        @range (0, 1)
    */
    typedef int boolean;

    /* An X/Y/Z style reference
    */
    typedef string obj_ref;

    /*
        required params:
        assembly_ref: Genome assembly object reference
        binned_contig_name: BinnedContig object name and output file header
        workspace_name: the name of the workspace it gets saved to.
        reads_list: list of reads object (PairedEndLibrary/SingleEndLibrary) upon which CONCOCT will be run

        optional params:
        thread: number of threads; default 1
        min_contig_length: minimum contig length; default 1000bp
        contig_split_size: length to split long contigs; default 10000bp
        ref: https://github.com/BinPro/CONCOCT

    */
    typedef structure {
        obj_ref assembly_ref;
        string binned_contig_name;
        string workspace_name;
        list<obj_ref> reads_list;

        int thread;
        int min_contig_length;
        int contig_split_size;
    } ConcoctInputParams;

    /*
        result_folder: folder path that holds all files generated by run_kb_concoct
        report_name: report name generated by KBaseReport
        report_ref: report reference generated by KBaseReport
    */
    typedef structure{
        string result_directory;
        obj_ref binned_contig_obj_ref;
        string report_name;
        string report_ref;
    }ConcoctResult;

    funcdef run_kb_concoct(ConcoctInputParams params)
        returns (ConcoctResult returnVal) authentication required;

};
