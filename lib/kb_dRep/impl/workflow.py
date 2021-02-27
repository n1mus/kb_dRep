

from .params import Params



def do_workflow(params):

    params = Params(params)

    app.update({
        'run_dir': os.path.join(app.shared_folder, 'run_kb_dRep_' + str(uuid.uuid4())), # folder dedicated to this API-method run
        'warnings': [],
        'params': params,
    })

    os.mkdir(app.run_dir)



    #
    ##
    ### load files, obj
    ####
    #####

    objs = []
    for upa in params['obj_upas']:

        dtype = Var.ws.get_object_info3({'objects': [{'ref': params['input_upa']}]})['infos'][0][2]

        if dtype.startswith('KBaseMetagenomes.BinnedContigs'):
            obj = BinnedContigs(upa)

        elif dtype.startswith('KBaseSets.GenomeSet'):
            obj = GenomeSet(upa)
            
        elif dtype.startswith('KBaseSets.AssemblySet'):
            obj = AssemblySet(upa)
            
        elif dtype.startswith('KBaseGenomes.Genome'):
            obj = Genome(upa)

        elif dtype.startswith('KBaseGenomeAnnotations.Assembly'):
            obj = Assembly(upa)

        else:
            raise Exception()

        objs.append(obj)



    #
    ##
    ### pool
    ####
    #####


    pooled_bins_dir = os.path.join(app.run_dir, 'pooled_bins')
    os.mkdir(pooled_bins_dir)

    for obj in obj:
        obj.pool_into(pooled_bins_dir)

    pooled_bins_fn_l = os.listdir(pooled_bins_dir)
    pooled_bins_fp_l = [os.path.join(pooled_bins_dir, fn) for fn in os.listdir(pooled_bins_dir)]

    dprint("os.listdir(pooled_bins_dir)")


    #
    ##
    ###
    #### params
    #####


    dRep_params_l = app.params.get_non_default_params_l()

    #
    if 'processors' in params:
        num_proc = params['processors'] # non-Narrative calls
    else: # TODO detect environment to assign different default for Narrative / non-Narrative
        num_proc = 8 # Narrative calls



    #
    ##
    ### run dRep
    ####
    #####


    dRep_work_dir = os.path.join(app.run_dir, 'dRep_work_dir')

    dRep_cmd = ([
        'dRep',
        'dereplicate',
        dRep_work_dir,
        '--genomes'
    ]
    + pooled_bins_fn_l
    + dRep_params_l
    + [
        '--debug',
        '--processors',
        str(num_proc),
    ])

    dRep_cmd = ' '.join(dRep_cmd)

    try:
        run_check(dRep_cmd)

    except Exception as e:

        # check: if Bdb.csv is empty -> nothing passed length/qual filtering
        with open(os.path.join(dRep_work_dir, 'data_tables/Bdb.csv')) as f:
            num_lines = sum(1 for line in f)

        if num_lines == 1:
            raise Exception('No assembly, bin, or genome passed filtering')

        else:
            raise(e)


    #app.dRep_cmd = dRep_cmd


    #
    ##
    ### save result set objs 
    ####
    #####

    derep_l = os.listdir(
        os.path.join(dRep_work_dir, 'dereplicated_genomes')
    )
    objects_created = []

    if params.getd('output_as_assembly'):
        assembly_ref_l = []
        for obj in objs:
            assembly_ref_l.append(
                obj.get_derep_assembly_refs()
            )
        ref = AssemblySet(ref_l=assembly_ref_l).save()

        objects_created.append(
            dict(
                ref=ref,
                description='Dereplication results',
            )
        )

    else:
        genome_l, assembly_l, genome_set_l, assembly_set_l, binned_contigs_l = partition_by_type(objs)

        # genome types
        genome_ref_l = aggregate_member_refs(genome_l, genome_set_l)
        ref = GenomeSet(ref_l=genome_ref_l).save()
        objects_created.append(
            dict(
                ref=ref,
                description='Dereplication results',
            )
        )

        # assembly types
        assembly_ref_l = aggregate_member_refs(assembly_l, assembly_set_l)
        ref = AssemblySet(ref_l=assembly_ref_l).save()
        objects_created.append(
            dict(
                ref=ref,
                description='Dereplication results',
            )
        )

        # binned contigs types
        for binned_contigs in binned_contigs_l:
            ref = binned_contigs.save_dereplicated()
            objects_created.append(
                dict(
                    ref=ref
                    description='Dereplication results',
                )
            )





    #
    ##
    ### html
    ####
    #####

    hb = report.HTMLBuilder(BinnedContigs.created_instances, dRep_cmd_str, dRep_work_dir)
    hb.build()
    html_dir, html_fp = hb.write()




    #
    ##
    ###
    ####
    #####

    file_links = [{
        'path': dRep_work_dir,
        'name': 'dRep_work_directory.zip',
        'description': 'Results (genomes, figures, etc.), intermediate files, logs, warnings ...',
    }]

    html_links = [{
        'path': html_dir,
        'name': os.path.basename(html_fp),
    }]

    report_params = {
            'warnings': app.warnings,
            'direct_html_link_index': 0,
            'html_links': html_links,
            'file_links': file_links,
            'report_object_name': 'kb_dRep_report',
            'workspace_id': params['workspace_id'],
            'objects_created': objects_created
            }

    if app.debug and params.get('skip_kbReport'):
        return

    report_output = app.kbr.create_extended_report(report_params)

    output = {
        'report_name': report_output['name'],
        'report_ref': report_output['ref'],
    }


    return output


def aggregate_derep_refs(star_l, star_set_l):
    '''
    :params star_l: genome or assembly list
    :params star_set_l: genome set or assembly set list
    '''

    star_ref_l = [
        star.ref for star in star_l
    ] + [
        ref for ref in
        star_set.derep_member_ref_l
        for star_set in star_set_l
    ]

    return star_ref_l

def partition_by_type(objs):
    genome_l = [obj for obj in objs if obj.TYPE == 'KBaseGenomes.Genome']
    assembly_l = [obj for obj in objs if obj.TYPE == 'KBaseGenomeAnnotations.Assembly']
    genome_set_l = [obj for obj in objs if obj.TYPE == 'KBaseSets.GenomeSet']
    assembly_set_l = [obj for obj in objs if obj.TYPE == 'KBaseSets.AssemblySet']
    binned_contigs_l = [obj for obj in objs if obj.TYPE == 'KBaseMetagenomes.BinnedContigs']

    return genome_l, assembly_l, genome_set_l, assembly_set_l, binned_contigs_l
