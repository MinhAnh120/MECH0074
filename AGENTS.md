# MAN - MECH0074

## Scope and context
- Use only `MinhAnh120/MECH0074` for all GitHub work. Develop on `stage1-pressure-waves`; do not merge without a request.
- This is MAN's Part B Stage 1 work, UCL MECH0074 Engineering in Extreme Environments, 2026/27.
- Read `docs/PROJECT_CONTEXT.md` and `docs/SOURCES.md` before changing the physical model. The brief controls assessment requirements; lecture notes and the supplied notebook control teaching conventions. Flag disagreements instead of silently changing the specification.
- Persist new decisions, unresolved questions and executed verification in the repository. Do not assume access to earlier chats or uncommitted attachments.

## Physics and implementation
- Adapt the supplied head/volume-flow MOC formulation with quasi-steady Darcy friction and a consistent steady initial head gradient.
- Use SI internally. Distinguish absolute pressure, pressure above external pressure, and change from the initial local pressure.
- Preserve prescribed linear velocity/flow closure over 0.020 s. Valve-area closure is a different model.
- Distinguish the supplied bulk sound speed (1100 m/s) from the lecture's effective pipe wave speed. Keep both labelled reference cases until the specification interpretation is settled.
- Darcy factor 0.020 remains a teaching comparison. Read docs/FRICTION_AND_PHASE.md for the provisional Colebrook nominal/bounds and pure-methane phase screen; these are not measured LNG or pipe properties. Preserve the distinction between published correlations, adopted assumptions and actual unknown composition/finish.
- Retain signed `Q*abs(Q)` friction, simultaneous old-to-new updates and `dt=dx/a`. Document any station interpolation.
- A mathematical single-phase verification trace is not an approved LNG design load. Report negative absolute pressure and missing/entered phase limits. Never clip pressure to conceal failure.
- Actual loop geometry and bend station coordinates are not supplied. Do not invent them, or label benchmark stations as bend histories.

## Working rules
- Use small, readable Python functions and restrained comments explaining physics or non-obvious choices. NumPy, Matplotlib and standard-library unittest suffice for this stage.
- Verify with independent analytical waves, steady flow, boundary conditions, the worked arithmetic, FFT frequencies and resolution comparisons. Run `python -m unittest discover -s tests -v` and `python complete_stage1.py --out results/latest` after relevant changes. Preserve legacy numerical verification as appropriate.
- Record actual pass/fail results, dependency versions and limitations. Do not claim experimental or code-compliance validation from numerical agreement.
- Make sequential, meaningful commits. Preserve other contributors' work and inspect the current branch head before updating it.
- When delivering produced or requested code, put each complete code/configuration file directly in a fenced copyable field in the chat. Give its repository path and exact commit. Links supplement these fields.
- Attribute supplied versus adapted versus new work. External scientific claims require peer-reviewed literature; identify software documentation and standards by their actual source type.
- Keep original course PDFs, archives and the full teaching notebook out of this public repository unless redistribution is specifically authorised. Cite their filenames, pages/cells and hashes in `docs/SOURCES.md`.
