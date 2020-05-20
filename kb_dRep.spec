/* A KBase module: kb_dRep
 */

module kb_dRep {

    /* Workspace reference in the form D/D/D
     * @id ws
     */
    typedef string ws_ref;

    /* 'True' or 'False'
     */
    typedef string bool;

    /* All optional
     * (same as original dRep filtering parameters)     
     */    
    typedef structure {
        int length;
        float completeness;
        float contamination;
        bool ignoreGenomeQuality;
    } params_filtering;

    /* All optional
     * (same as original dRep genome_comparison parameters)     
     */    
    typedef structure {
        int MASH_sketch;
        string S_algorithm;
        string n_PRESET;
    } params_genome_comparison;   

    /* All optional
     * (same as original dRep clustering parameters)     
     */    
    typedef structure {
        float P_ani;
        float S_ani;
        bool SkipMash;
        bool SkipSecondary;
        float cov_thresh;
        string coverage_method;
        string clusterAlg;
    } params_clustering;

    /* All optional
     * (same as original dRep scoring parameters)     
     */    
    typedef structure {
        float completeness_weight;
        float contamination_weight;
        float strain_heterogeneity_weight;
        float N50_weight;
        float size_weight;
    } params_scoring;

    /* All optional
     * (same as original dRep warnings parameters)     
     */    
    typedef structure {
        float warn_dist;
        float warn_sim;
        float warn_aln;
    } params_warnings;
    
    /* All optional except `genomes_refs` and `workspace_id`
     * (mostly same parameters and parameter groupings as original dRep parameters)     
     * 
     * Parameters are grouped here in accordance with how the narrative passes the parameters
     * However, parameters are also accepted in flattened format, e.g., during direct app calls
     *
     * If a mixture of flattened and grouped parameters is used, no behavior is guaranteed
     */
    typedef structure {
        list<ws_ref> genomes_refs;
        params_filtering filtering;
        params_genome_comparison genome_comparison;
        params_clustering clustering;
        params_scoring scoring;
        params_warnings warnings;
        string checkM_method;
        string workspace_name;
        int workspace_id;
    } params_dereplicate;

    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /* Given BinnedContigs references and parameters, pool and dereplicate the bins
     * Create new BinnedContigs objects named `original_name + '.dRep'`
     * Time and memory (40GB) intensive
     */   
    funcdef run_dereplicate(params_dereplicate params) 
        returns (ReportResults output) authentication required;

};
