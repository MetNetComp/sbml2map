# Escher viewer and validation files

- `map-schema-1-0-0.json`: Official Escher map JSON schema.
- `escher-1.8.1.min.js`: Escher 1.8.1 viewer, embedded in offline previews.
- `LICENSE`: MIT license for the bundled library.
- `sources.json`: Source URLs and SHA-256 hashes.

The schema comes from the [official repository](https://github.com/zakandrewking/escher/blob/master/jsonschema/1-0-0).
The viewer comes from the [versioned npm package](https://unpkg.com/escher@1.8.1/dist/escher.min.js).
Maps use the Escher 1-0-0 format. Output reports record the viewer version used for browser checks.
