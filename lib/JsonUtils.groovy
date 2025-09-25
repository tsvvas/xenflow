class JsonUtils {
  static String toPrettyJson(def obj) {
    return groovy.json.JsonOutput.prettyPrint(
      groovy.json.JsonOutput.toJson(obj)
    )
  }
}
