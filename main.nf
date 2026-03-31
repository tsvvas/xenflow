nextflow.enable.dsl = 2

include { CONVERT_XENIUM }      from './modules/convert_xenium.nf'
include { RESEGMENT_NUCLEI }    from './modules/resegment_nuclei.nf'
include { RESEGMENT_CELLS }     from './modules/resegment_cells.nf'
include { DETECT_TISSUE }       from './modules/detect_tissue.nf'
include { SPLIT_SAMPLES }       from './modules/split_samples.nf'
include { MAP_REFERENCE }       from './modules/map_reference.nf'
include { IDENTIFY_PROGRAMS }   from './modules/identify_programs.nf'

workflow {
    channel
        .fromPath( params.xenium_dir, checkIfExists: true )
        .map { p -> tuple( p, p.name.split('__')[1] ) }
        .set { xenium_ch }

    convert_out       = CONVERT_XENIUM( xenium_ch )
    nuclei_reseg_out  = RESEGMENT_NUCLEI( convert_out )
    if( params.resegment_cells ) {
            cells_reseg_out = RESEGMENT_CELLS( nuclei_reseg_out )
    }
        else {
            cells_reseg_out = nuclei_reseg_out
    }
    detect_out        = DETECT_TISSUE( cells_reseg_out )

    regions_ch     = detect_out.map { _bbox, regions, _fig, id -> tuple( regions, id ) }
    convert_keyed  = cells_reseg_out.map { z, id -> tuple( id, z ) }
    regions_keyed  = regions_ch.map { r, id -> tuple( id, r ) }

    split_in   = convert_keyed.join( regions_keyed ).map { id, z, r -> tuple( z, r, id ) }
    split_out  = SPLIT_SAMPLES( split_in )

    map_in = split_out.flatMap { h5ads, id ->
        h5ads.collect { h5ad -> tuple(h5ad, id) }
    }
    _map_out = MAP_REFERENCE( map_in )

    _identify_out = IDENTIFY_PROGRAMS( map_in )
}


workflow TEST {
    xenium_ch = channel
        .fromPath( params.test_zarr, checkIfExists: true )
        .map { zarr ->
            def name = zarr.baseName
            def idx = name.lastIndexOf('_')
            def sample_id = name.substring(0, idx)
            tuple( zarr, sample_id )
        }

    convert_out      = CONVERT_XENIUM( xenium_ch )
    nuclei_reseg_out = RESEGMENT_NUCLEI( convert_out )
    if( params.resegment_cells ) {
            cells_reseg_out = RESEGMENT_CELLS( nuclei_reseg_out )
    }
        else {
            cells_reseg_out = nuclei_reseg_out
    }
    detect_out        = DETECT_TISSUE( cells_reseg_out )
    regions_ch     = detect_out.map { _bbox, regions, _fig, id -> tuple( regions, id ) }
    convert_keyed  = cells_reseg_out.map { z, id -> tuple( id, z ) }
    regions_keyed  = regions_ch.map { r, id -> tuple( id, r ) }

    split_in   = convert_keyed.join( regions_keyed ).map { id, z, r -> tuple( z, r, id ) }
    split_out  = SPLIT_SAMPLES( split_in )

    map_in = split_out.flatMap { h5ads, id ->
        h5ads.collect { h5ad -> tuple(h5ad, id) }
    }
    _map_out = MAP_REFERENCE( map_in )

    _identify_out = IDENTIFY_PROGRAMS( map_in )
}
