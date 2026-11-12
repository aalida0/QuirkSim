# Bugfix: MG5 3.5.1 missing parameter-card generator

MG5 imported the full model and constructed 34 partonic processes, but event
generation stopped before integration with:

```
can't open .../bin/internal/ufomodel/write_param_card.py
FileNotFoundError: .../bin/internal/ufomodel/param_card.dat
```

The self-contained model previously omitted `write_param_card.py`, which
MG5 executes at the start of a run. This patch vendors that file from the
upstream MadGraph5_aMC Standard Model UFO `models/sm/write_param_card.py`
(3.x branch) and changes our builder to retain the writer.

The corrected MG5 script no longer tries to set an unavailable shower
switch and sets the charged mass once, through the shared `MQC` parameter.

After updating the branch, rerun the `examples/pp_dirac_triplet_13p6TeV.mg5`
script. Its `output ... -f` will overwrite a same-name output directory,
so preserve useful existing results before re-running.

Status: event-generation validation is tracked separately in GitHub Actions.
