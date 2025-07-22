nextflow.enable.dsl = 2

include { CONVERT_XENIUM }     from './modules/convert_xenium.nf'
include { DETECT_TISSUE }      from './modules/detect_tissue.nf'
include { SPLIT_SAMPLES }      from './modules/split_samples.nf'
include { IDENTIFY_PROGRAMS }  from './modules/identify_programs.nf'

workflow {
    Channel
        .fromPath( params.xenium_dir, checkIfExists: true )
        .map { p -> tuple( p, p.name.split('__')[1] ) } 
        .set { xenium_ch }
    convert_out = CONVERT_XENIUM( xenium_ch )
    detect_out  = DETECT_TISSUE( convert_out )

    regions_ch = detect_out.map { bbox, regions, fig, id -> tuple( regions, id ) }
    convert_keyed = convert_out.map { z, id -> tuple( id, z ) }
    regions_keyed = regions_ch .map { r, id -> tuple( id, r ) }

    split_in = convert_keyed.join( regions_keyed ).map { id, z, r -> tuple( z, r, id ) }
    split_out = SPLIT_SAMPLES( split_in )

    identify_in = split_out.flatMap { path, _ -> path }
    identify_out = IDENTIFY_PROGRAMS( identify_in )
}
