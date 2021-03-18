/* A KBase module: kb_dRep
 */

module kb_dRep {

    /*
     * @id ws
     */
    typedef string ws_ref;

    /* 0 or 1
     */
    typedef int bool;

    typedef structure {
        int length;
        float completeness;
        float contamination;
    } filtering;

    typedef structure {
        int MASH_sketch;
        string S_algorithm;
        string n_PRESET;
    } genome_comparison;   

    typedef structure {
        float P_ani;
        float S_ani;
        float cov_thresh;
        string coverage_method;
        string clusterAlg;
    } clustering;

    typedef structure {
        float completeness_weight;
        float contamination_weight;
        float strain_heterogeneity_weight;
        float N50_weight;
        float size_weight;
    } scoring;

    /* All optional except `obj_refs`, `workspace_name`, and `workspace_id`
     */
    typedef structure {
        list<ws_ref> obj_refs;
        filtering filtering;
        genome_comparison genome_comparison;
        clustering clustering;
        scoring scoring;
        string checkM_method;
        int processors;
        bool output_as_assembly;
        string output_suffix;
        string workspace_name;
        int workspace_id;
    } dRepParams;

    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    funcdef run_dereplicate(dRepParams params) 
        returns (ReportResults output) authentication required;

};
