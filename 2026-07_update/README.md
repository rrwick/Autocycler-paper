Two new versions and one new assembler to test:
- [metaMDBG](https://github.com/GaetanBenoitDev/metaMDBG) v1.4
- [Myloasm](https://github.com/bluenote-1577/myloasm) v0.6.0
- [Hybracter](https://github.com/gbouras13/hybracter) v0.14.0
- [Ilesta](https://github.com/yvlaere/Ilesta) v1.2.1

Not running these assemblers/piplines, since their version hasn't changed since the Autocycler paper:
- [Canu](https://github.com/marbl/canu) 2.3
- [Flye](https://github.com/mikolmogorov/Flye) 2.9.6
- [Hifiasm](https://github.com/chhylp123/hifiasm) v0.25.0
- [LJA](https://github.com/AntonBankevich/LJA/) v0.2
- [Miniasm](https://github.com/lh3/miniasm) 0.3
  * [Minipolish](https://github.com/rrwick/Minipolish) has had some updates (v0.1.3 to 0.2.1), but nothing significant enough to warrant retesting.
- [NECAT](https://github.com/xiaochuanle/NECAT) 0.0.1_update20200803
- [NextDenovo](https://github.com/nextomics/NextDenovo) 2.5.2, [NextPolish](https://github.com/Nextomics/NextPolish) 1.4.1
- [Raven](https://github.com/lbcb-sci/raven) 1.8.3
- [Wtdbg](https://github.com/ruanjue/wtdbg2) 2.5
- [Dragonflye](https://github.com/rpetit3/dragonflye) v1.2.1

I'll then try two Autocycler pipelines:
- [`autocycler_full.sh`](https://github.com/rrwick/Autocycler/tree/main/pipelines/Automated_Autocycler_Bash_script_by_Ryan_Wick): this one was run for the Autocycler paper, but it's had a few little updates (e.g. the inclusion of Myloasm).
- [`autocycler_full_fast.sh`](): this is a new one with some tweaks to cut down execution time.
