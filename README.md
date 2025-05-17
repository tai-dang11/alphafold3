# Protenix: Protein + X

A trainable PyTorch reproduction of [AlphaFold 3](https://www.nature.com/articles/s41586-024-07487-w).

![Protenix predictions](assets/protenix_predictions.gif)

## 🔥 Feature Update
* 🚀 The preview version of [constraint feature](./README.md#early-access-to-new-constraint-feature) is released to branch [`constraint_esm`](https://github.com/bytedance/Protenix/tree/constraint_esm).
* 🪐 The [training data pipeline](./docs/prepare_training_data.md) is released.
* ⚡️  The [MSA pipeline](./docs/msa_pipeline.md) is released.
* 🛸 Use [local colabfold_search](./docs/colabfold_compatiable_msa.md) to generate protenix-compatible MSA.

## Installation

### Run with PyPI (recommended):

```bash
pip3 install protenix
```

## Inference

### Command line inference

If you set up `Protenix` by `pip`, you can run the following command to do model inference:

```bash
# the default n_cycle/n_step/n_samples is 10/200/5 respectively, you can modify it by passing --cycle x1 --step x2 --sample x3

# run with example.json, which contains precomputed msa dir.
protenix predict --input examples/example.json --out_dir  ./output --seeds 101

# run with multiple json files, the default seed is 101.
protenix predict --input ./jsons_dir/ --out_dir  ./output

# if the json do not contain precomputed msa dir,
# add --use_msa_server to search msa and then predict.
# if mutiple seeds are provided, split them by comma.
protenix predict --input examples/example_without_msa.json --out_dir ./output --seeds 101,102 --use_msa_server
```

### Convert PDB/CIF file to json

If your input is pdb or cif file, you can convert it to json file for inference.

```bash
# ensure `release_data/ccd_cache/components.cif` or run:
python scripts/gen_ccd_cache.py -c release_data/ccd_cache/ -n [num_cpu]

# for PDB
# download pdb file
wget https://files.rcsb.org/download/7pzb.pdb
# run with pdb/cif file, and convert it to json file for inference.
protenix tojson --input examples/7pzb.pdb --out_dir ./output

# for CIF (same process)
# download cif file
wget https://files.rcsb.org/download/7pzb.cif
# run with pdb/cif file, and convert it to json file for inference.
protenix tojson --input examples/7pzb.cif --out_dir ./output
```

### Performance details

**Detailed information on the format of the input JSON file and the output files can be found in [<u> input and output documentation </u>](docs/infer_json_format.md)**.

Alternatively you can run inference by:

Note: by default, we do not use layernorm and EvoformerAttention kernels for simple configuration, if you want to speed up inference, see [<u> setting up kernels documentation </u>](docs/kernels.md).

```bash
bash inference_demo.sh
```

Arguments in this scripts are explained as follows:

* `input_json_path`: path to a JSON file that fully describes the input.
* `dump_dir`: path to a directory where the results of the inference will be saved.
* `dtype`: data type used in inference. Valid options include `"bf16"` and `"fp32"`.
* `use_msa`: whether to use the MSA feature, the default is true.
* `use_esm`: whether to use the ESM feature, the default is false.


### Convert PDB/CIF file to json

If your input is pdb or cif file, you can convert it to json file for inference.
```bash
# run with pdb/cif file, and convert it to json file for inference.
protenix tojson --input examples/7pzb.pdb --out_dir ./output
```

### MSA search
We also provide an independent MSA search function, you can do msa search from json file or fasta file.
```bash
# run msa search with json file, it will write precomputed msa dir info to a new json file.
protenix msa --input examples/example_without_msa.json --out_dir ./output

# run msa search with fasta file which only contains protein.
protenix msa --input examples/prot.fasta --out_dir ./output
```

## Training
If you're interested in model training, see [<u> training documentation </u>](docs/training.md).

## Performance
#### **Model Performance across Several Benchmarks**
![Overall Metrics](assets/overall_metrics.png)

#### ***Early Access to NEW Constraint Feature!***

🎉 Protenix now allows users to specify ***contacts***, enabling the model to leverage additional inter-chain information as constraint guidance! We benchmarked our constraint feature on Posebuster and a protein-antibody interfaces subset. Protenix demonstrates powerful ability in predicting more accurate structures. If you want to have a try, checkout to branch `constraint_esm` for details about the input format.

![Constraint Metrics](assets/constraint_metrics.png)

> **Tips:** Our online service already supports the new features, so feel free to try it now! Due to the preview version, the constraint support is only applicable in the branch `constraint_esm`. If you want to run inference via the command line, please check out to this branch first.

## Training and Inference Cost

See the [<u>model_train_inference_cost documentation</u>](docs/model_train_inference_cost.md) for memory and time consumption in training and inference.


## Citing This Work

If you use this code or the model in your research, please cite the following paper:

```
@article{chen2025protenix,
  title={Protenix - Advancing Structure Prediction Through a Comprehensive AlphaFold3 Reproduction},
  author={Chen, Xinshi and Zhang, Yuxuan and Lu, Chan and Ma, Wenzhi and Guan, Jiaqi and Gong, Chengyue and Yang, Jincai and Zhang, Hanyu and Zhang, Ke and Wu, Shenghao and Zhou, Kuangqi and Yang, Yanping and Liu, Zhenyu and Wang, Lan and Shi, Bo and Shi, Shaochen and Xiao, Wenzhi},
  year={2025},
  doi = {10.1101/2025.01.08.631967},
  journal = {bioRxiv}
}
```

## Acknowledgements

Implementation of the layernorm operators referred to [OneFlow](https://github.com/Oneflow-Inc/oneflow) and [FastFold](https://github.com/hpcaitech/FastFold). We used [OpenFold](https://github.com/aqlaboratory/openfold) for some [module](protenix/openfold_local/) implementations, except the [`LayerNorm`](protenix/model/layer_norm/).