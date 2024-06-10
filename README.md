[![Build Status](https://travis-ci.org/WormBase/genedesc_generator.svg?branch=master)](https://travis-ci.org/WormBase/genedesc_generator) [![Coverage Status](https://coveralls.io/repos/github/valearna/wb_genedescriptions/badge.svg?branch=master&service=github)](https://coveralls.io/github/valearna/wb_genedescriptions?branch=master) [![Codacy Badge](https://api.codacy.com/project/badge/Grade/7a999c3a60f44df9a0312fdab82e405c)](https://www.codacy.com/app/valearna/wb_genedescriptions?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=valearna/wb_genedescriptions&amp;utm_campaign=Badge_Grade) [![Documentation Status](https://readthedocs.org/projects/wb-genedescriptions/badge/?version=latest)](http://wb-genedescriptions.readthedocs.io/en/latest/?badge=latest) [![Maintainability](https://api.codeclimate.com/v1/badges/a444e953e3dceab54e25/maintainability)](https://codeclimate.com/github/valearna/genedesc_generator/maintainability)

## Introduction

This software generates gene summaries for Model Organisms part of the Alliance of Genome Resources project.

## Installation

To install the software, clone the repository and run the following command from the
root directory of the project:

```bash
$ pip3 install .
```

## Generating summaries for WormBase

To generate gene summaries for WormBase, run the `exec_all_pipelines.sh` shell script in the `wormbase` folder.

### Configuring the main script

the main script can be configured through a configuration file in .yaml format. The default file is *genedesc.yaml*,
but a different file can be provided by using the argument -c.

To configure the WormBase pipeline, modify the file `config_wb.yml` in the `wormbase` folder.
