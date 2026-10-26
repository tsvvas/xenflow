import java.nio.file.Files
import java.nio.file.Path

class ConfigValidator {
    static void validate(List<Path> xeniumBundles, Path metadataFile, Path referenceFile, String projectDir) {
        def errors = []

        if (!projectDir?.trim()) {
            errors << 'PROJECTDIR must be set to a nonempty value.'
        }

        xeniumBundles.each { bundle ->
            if (!Files.isRegularFile(bundle.resolve('experiment.xenium')) &&
                !Files.isRegularFile(bundle.resolve('outs/experiment.xenium'))) {
                errors << "Xenium bundle '${bundle}' must contain experiment.xenium or outs/experiment.xenium."
            }
        }

        if (metadataFile != null && !Files.isRegularFile(metadataFile)) {
            errors << "Custom metadata file does not exist or is not a file: '${metadataFile}'."
        }
        if (referenceFile != null && !Files.isRegularFile(referenceFile)) {
            errors << "Reference dataset does not exist or is not a file: '${referenceFile}'."
        }

        if (errors) {
            throw new IllegalArgumentException('Invalid configuration:\n - ' + errors.join('\n - '))
        }
    }
}