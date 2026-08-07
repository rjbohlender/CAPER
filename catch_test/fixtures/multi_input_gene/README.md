# Multi-input single-gene fixture

This fixture exercises the multi-input `--input` mode for gene `SPLIT1`.

- `split_1.matrix` and `split_2.matrix` contain different portions of `SPLIT1`.
- Unrelated genes interrupt `SPLIT1` in both input files.
- Transcript `SPLIT1.2` occurs before `SPLIT1.1` in the first input, requiring the collected rows to be sorted before constructing the gene.
- `split_2.matrix` uses sample order `control3, case2, control1, case4, control4, case1, control2, case3`.
- `merged.matrix` contains the equivalent target rows in canonical sample and transcript order.

From the repository root, compare the modes with:

```bash
mkdir -p /tmp/caper-multi /tmp/caper-merged

./build/caper \
  -i catch_test/fixtures/multi_input_gene/split_1.matrix \
     catch_test/fixtures/multi_input_gene/split_2.matrix \
  -l SPLIT1 \
  -p catch_test/fixtures/multi_input_gene/samples.ped \
  -o /tmp/caper-multi \
  -m BURDEN --no_weights --nperm 1000 --successes 1001 --seed 8675309 -t 2

./build/caper \
  -i catch_test/fixtures/multi_input_gene/merged.matrix \
  -l SPLIT1 \
  -p catch_test/fixtures/multi_input_gene/samples.ped \
  -o /tmp/caper-merged \
  -m BURDEN --no_weights --nperm 1000 --successes 1001 --seed 8675309 -t 2

# The first line records the command and therefore has different input/output paths.
diff -u \
  <(tail -n +2 /tmp/caper-multi/BURDEN.simple) \
  <(tail -n +2 /tmp/caper-merged/BURDEN.simple)
diff -u /tmp/caper-multi/BURDEN.detail /tmp/caper-merged/BURDEN.detail
```

Both `diff` commands should produce no output.
