import java.nio.file.Files
import nextflow.Nextflow

class ConfigValidator {
    static void validate(Map params) {
        def errors = []

        if (!System.getenv('PROJECTDIR')?.trim()) {
            errors << 'PROJECTDIR must be set to a nonempty value.'
        }

        if (params.xenium_dir) {
            def inputs = params.xenium_dir instanceof List ? params.xenium_dir : [params.xenium_dir]
            inputs.each { input ->
                Nextflow.files([type: 'dir', checkIfExists: true], input).each { bundle ->
                    if (!Files.isRegularFile(bundle.resolve('experiment.xenium')) &&
                        !Files.isRegularFile(bundle.resolve('outs/experiment.xenium'))) {
                        errors << "Xenium bundle '${bundle}' must contain experiment.xenium or outs/experiment.xenium."
                    }
                }
            }
        }

        if (params.map_reference && !params.reference_dataset) {
            errors << 'Reference mapping requires --reference_dataset.'
        }

        [metadata_file: 'Custom metadata file', reference_dataset: 'Reference dataset'].each { key, label ->
            def value = params[key]
            if (value && !Files.isRegularFile(Nextflow.file([glob: false], value))) {
                errors << "${label} does not exist or is not a file: '${value}'."
            }
        }

        if (errors) {
            throw new IllegalArgumentException('Invalid configuration:\n - ' + errors.join('\n - '))
        }
    }
}
