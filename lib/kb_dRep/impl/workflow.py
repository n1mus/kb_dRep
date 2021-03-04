import os
import logging
import uuid

from .params import Params
from .config import app
from .kb_obj import Assembly, AssemblySet, Genome, GenomeSet, BinnedContigs
from . import report
from ..util.cli import run_check
from ..util.debug import dprint



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
    ### directories
    ####
    #####

    dRep_dir = os.path.join(app.run_dir, 'dRep_dir')
    report_dir = os.path.join(app.run_dir, 'report')

    app.update(
        dict(
            dRep_dir=dRep_dir,
            report_dir=report_dir,
        )
    )


    #
    ##
    ### load files, obj
    ####
    #####

    objs = []
    for upa in params['obj_refs']:

        type_ = app.ws.get_object_info3({
            'objects': [{'ref': upa}]
        })['infos'][0][2]

        if type_.startswith('KBaseMetagenomes.BinnedContigs'):
            obj = BinnedContigs(upa)

        elif type_.startswith('KBaseSets.GenomeSet'):
            obj = GenomeSet(upa, type='KBaseSets.GenomeSet')
 
        elif type_.startswith('KBaseSearch.GenomeSet'):
            obj = GenomeSet(upa, type='KBaseSearch.GenomeSet')
                    
        elif type_.startswith('KBaseSets.AssemblySet'):
            obj = AssemblySet(upa)
            
        elif type_.startswith('KBaseGenomes.Genome'):
            obj = Genome(upa)

        elif type_.startswith('KBaseGenomeAnnotations.Assembly'):
            obj = Assembly(upa)

        else:
            raise Exception(type_)

        objs.append(obj)



    #
    ##
    ### pool
    ####
    #####


    pooled_bins_dir = os.path.join(app.run_dir, 'pooled_bins')
    os.mkdir(pooled_bins_dir)

    for obj in objs:
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



    dRep_cmd = ([
        'dRep',
        'dereplicate',
        dRep_dir,
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
        with open(os.path.join(dRep_dir, 'data_tables/Bdb.csv')) as f:
            num_lines = sum(1 for line in f)

        if num_lines == 1:
            raise Exception('No assembly, bin, or genome passed filtering')

        else:
            raise(e)


    #
    ##
    ### save result set objs 
    ####
    #####

    derep_l = os.listdir(
        os.path.join(dRep_dir, 'dereplicated_genomes')
    )
    for obj in objs:
        obj.identify_dereplicated(derep_l)

    objects_created = []
    description = 'Dereplication results'

    if params.getd('output_as_assembly'):
        assembly_ref_l = []
        for obj in objs:
            assembly_ref_l.extend(
                obj.get_derep_assembly_refs()
            )
        assembly_ref_l = uniq(assembly_ref_l)
        ref = AssemblySet(ref_l=assembly_ref_l).save('Assemblies' + params.getd('output_suffix'))

        objects_created.append(
            dict(
                ref=ref,
                description=description,
            )
        )

    else:
        genome_l, assembly_l, genome_set_l, assembly_set_l, binned_contigs_l = partition_by_type(objs)

        # genome types
        genome_ref_l = aggregate_derep_element_refs(genome_l, genome_set_l)
        if len(genome_ref_l) > 0:
            ref = GenomeSet(ref_l=genome_ref_l).save('Genomes' + params.getd('output_suffix'))
            objects_created.append(
                dict(
                    ref=ref,
                    description=description,
                )
            )

        # assembly types
        assembly_ref_l = aggregate_derep_element_refs(assembly_l, assembly_set_l)
        if len(assembly_ref_l) > 0:
            ref = AssemblySet(ref_l=assembly_ref_l).save('Assemblies' + params.getd('output_suffix'))
            objects_created.append(
                dict(
                    ref=ref,
                    description=description,
                )
            )

        # binned contigs types
        for binned_contigs in binned_contigs_l:
            if not binned_contigs.is_fully_dereplicated():
                ref = binned_contigs.save_dereplicated(binned_contigs.name + params.getd('output_suffix'))
                objects_created.append(
                    dict(
                        ref=ref,
                        description=description,
                    )
                )





    #
    ##
    ### html
    ####
    #####

    hb = report.HTMLBuilder(dRep_dir, report_dir)
    report_fp = hb.write()




    #
    ##
    ###
    ####
    #####

    file_links = [{
        'path': dRep_dir,
        'name': 'results.zip',
    }]

    html_links = [{
        'path': report_dir,
        'name': os.path.basename(report_fp),
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

    report_output = app.kbr.create_extended_report(report_params)

    output = {
        'report_name': report_output['name'],
        'report_ref': report_output['ref'],
    }


    return output


def aggregate_derep_element_refs(star_l, star_set_l):
    '''
    :params star_l: genome or assembly list
    :params star_set_l: genome set or assembly set list
    '''

    star_ref_l = [
        star.ref for star in star_l
    ] + [
        ref 
        for star_set in star_set_l
        for ref in star_set.get_derep_member_refs() 
    ]

    return uniq(star_ref_l)

def partition_by_type(objs):
    genome_l = [obj for obj in objs if obj.TYPE == 'KBaseGenomes.Genome']
    assembly_l = [obj for obj in objs if obj.TYPE == 'KBaseGenomeAnnotations.Assembly']
    genome_set_l = [obj for obj in objs if obj.TYPE == 'KBaseSets.GenomeSet' or obj.TYPE == 'KBaseSearch.GenomeSet']
    assembly_set_l = [obj for obj in objs if obj.TYPE == 'KBaseSets.AssemblySet']
    binned_contigs_l = [obj for obj in objs if obj.TYPE == 'KBaseMetagenomes.BinnedContigs']

    return genome_l, assembly_l, genome_set_l, assembly_set_l, binned_contigs_l


def uniq(l):
    return sorted(list(set(l)))

