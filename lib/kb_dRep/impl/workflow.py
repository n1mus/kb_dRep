import os
import logging
import uuid

from .params import Params
from .config import app, ref_leaf
from .kb_obj import Assembly, AssemblySet, Genome, GenomeSet, BinnedContigs
from . import report
from ..util.cli import run_check
from ..util.debug import dprint



def do_workflow(params):

    params = Params(params)

    app.update({
        'run_dir': os.path.join(app.shared_folder, 'run_kb_dRep_' + str(uuid.uuid4())), # folder dedicated to this API-method run
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

        if type_.startswith(BinnedContigs.TYPE):
            obj = BinnedContigs(upa)

        elif type_.startswith(GenomeSet.TYPE) or type_.startswith(GenomeSet.LEGACY_TYPE):
            obj = GenomeSet(ref=upa)
 
        elif type_.startswith(AssemblySet.TYPE):
            obj = AssemblySet(ref=upa)
            
        elif type_.startswith(Genome.TYPE):
            obj = Genome(upa)

        elif type_.startswith(Assembly.TYPE):
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
    ### run dRep
    ####
    #####

    dRep_params_l = params.get_non_default_tool_params()

    dRep_cmd = ([
        'dRep dereplicate',
        dRep_dir,
        '--genomes'
    ]
    + pooled_bins_fp_l
    + dRep_params_l
    + [
        '--debug',
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

    objects_created = save_results(objs, params, dRep_dir)



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
        'direct_html_link_index': 0,
        'html_links': html_links,
        'file_links': file_links,
        'report_object_name': 'kb_dRep_report',
        'workspace_id': params['workspace_id'],
        'objects_created': objects_created,
    }

    report_output = app.kbr.create_extended_report(report_params)

    output = {
        'report_name': report_output['name'],
        'report_ref': report_output['ref'],
    }


    return output


def save_results(objs, params, dRep_dir):
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
            if obj.TYPE == BinnedContigs.TYPE:
                obj.save_derep_as_assemblies(params['workspace_name'])
            assembly_ref_l.extend(obj.get_derep_assembly_refs())
        assembly_ref_l = uniq_refs(assembly_ref_l)
        ref = AssemblySet(ref_l=assembly_ref_l).save(
            'Assemblies' + params.getd('output_suffix'),
            params['workspace_id']
        )

        objects_created.append(
            dict(
                ref=ref,
                description=description,
            )
        )

    else:
        assembly_l, genome_l, assembly_set_l, genome_set_l, binned_contigs_l = partition_by_type(objs)

        # assembly types
        assembly_ref_l = aggregate_derep_element_refs(assembly_l, assembly_set_l)
        if len(assembly_ref_l) > 0:
            ref = AssemblySet(ref_l=assembly_ref_l).save(
                'Assemblies' + params.getd('output_suffix'),
                params['workspace_id']
            )
            objects_created.append(
                dict(
                    ref=ref,
                    description=description,
                )
            )

        # genome types
        genome_ref_l = aggregate_derep_element_refs(genome_l, genome_set_l)
        if len(genome_ref_l) > 0:
            ref = GenomeSet(ref_l=genome_ref_l).save(
                'Genomes' + params.getd('output_suffix'),
                params['workspace_id']
            )
            objects_created.append(
                dict(
                    ref=ref,
                    description=description,
                )
            )

        # binned contigs types
        for binned_contigs in binned_contigs_l:
            if not binned_contigs.is_fully_dereplicated():
                ref = binned_contigs.save_dereplicated(
                    binned_contigs.name + params.getd('output_suffix'),
                    params['workspace_id']
                )
                objects_created.append(
                    dict(
                        ref=ref,
                        description=description,
                    )
                )

    return objects_created


def aggregate_derep_element_refs(star_l, star_set_l):
    """
    :params star_l: genome or assembly list
    :params star_set_l: genome set or assembly set list
    """

    star_ref_l = [
        star.ref for star in star_l if star.in_derep
    ] + [
        ref 
        for star_set in star_set_l
        for ref in star_set.get_derep_member_refs() 
    ]

    return uniq_refs(star_ref_l)


def partition_by_type(objs):
    assembly_l = [obj for obj in objs if obj.TYPE == Assembly.TYPE]
    genome_l = [obj for obj in objs if obj.TYPE == Genome.TYPE]
    assembly_set_l = [obj for obj in objs if obj.TYPE == AssemblySet.TYPE]
    genome_set_l = [obj for obj in objs if obj.TYPE == GenomeSet.TYPE or obj.TYPE == GenomeSet.LEGACY_TYPE]
    binned_contigs_l = [obj for obj in objs if obj.TYPE == BinnedContigs.TYPE]

    return assembly_l, genome_l, assembly_set_l, genome_set_l, binned_contigs_l


def uniq_refs(l):
    d = {} # ref leaf to shortest ref path

    for e in l:
        lf = ref_leaf(e)
        if lf not in d or e.count(';') < d[lf].count(';'):
            d[lf] = e 

    return list(d.values())

