# Parameter and other files

## Parameter file

The parameter file(s) are created in the [YAML](#yaml-files) format. All the
defaults are located in
the [`cg/default.yml`](https://github.com/vitreusx/pas-cg/tree/main/cg/src/default.yml)
master file.

## YAML files

The [YAML](https://en.wikipedia.org/wiki/YAML) files used throughout the program
are standard YAML files, with few "conventions" added. In the documentation, we
specify schemas for the entries, i.e. associate types with various nodes. Aside
from the usual types, like numbers or strings or lists etc., we introduce a
number of custom types.

### File nodes

Whenever a given YAML entry in the documentation is denoted as `file`, one can
either provide it as a scalar with the contents, or refer to it by a path. The
path is relative **to the YAML file**.

```yaml
file: scalar

file:
  (at path): relpath
```

This is chiefly done in order to not have to include the file verbatim in the
YAML file, while ensuring the consistency of such "external" references.

### Subnodes

In certain cases an entry is denoted as `subnode` - this means that either we
can provide its value as a regular YAML subnode, or refer to it in another file
in a manner as above.

### CSV files

A CSV file (usually denoted as `csv file`) can be included in the YAML file in
the same way that a regular file is, by contents or by a path.