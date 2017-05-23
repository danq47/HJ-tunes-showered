<TeXmacs|1.99.2>

<style|generic>

<\body>
  <doc-data|<doc-title|Mass effects>|<doc-author|<author-data|<author-name|K.
  Hamilton, P. Nason and G. Zanderighi>>>>

  This manual documents the upgrade of the <verbatim|HJ-MiNLO> code
  illustrated in ref. <cite|masseffects> <cite|Hamilton:2013fea>
  <cite|Hamilton:2012rf>, where we describe the implementation of quark-mass
  effects in the <with|font-family|tt|HJ-MiNLO> generator, and its NNLO
  accurate reweighting. It also updates the setup found in
  <verbatim|HJ/NNLOPS>, describings the reweighting of the
  <verbatim|HJ-MiNLO> ouptut in order to achieve NNLO accuracy, since in the
  new code we use extensively the reweighting conventions agreed upon in the
  Les Houches convention described in <hlink|http://phystev.in2p3.fr/wiki/2013:groups:tools_lheextension|>.
  The new NNLOPS code is now in <verbatim|HJ/NNLOPS-mass-effects> and in the
  main <verbatim|HJ> directory. The <verbatim|HJ/NNLOPS> directory has not
  been modified, and the manual included there is still useful for all
  aspects that have not been updated in the present manual.\ 

  The use of the code is illustrated in the directory
  <with|font-family|tt|HJ/PaperRun>, where a script for reproducing the runs
  performed for the preparation of ref. <cite|masseffects> is given. Here we
  refer to this example in order to illustrate how the program works.

  <section|Production of the HJ-MiNLO events with mass effects.>

  \;

  The inclusion of mass effects is controlled by the following flags in the
  <with|font-family|tt|powheg.input> file:

  <\verbatim-code>
    quarkmasseffects 1 \ \ ! Include quark mass effects

    topmass 172.5d0 \ \ \ \ \ ! top quark mass (172.5 in
    HNNLO-patches/mdata.f)

    bottommass 3.38d0 \ \ \ ! bottom quark mass in MSbar at MH

    bmass_in_minlo 1 \ \ \ \ ! Include quark mass effects in Sudakov exponent
  </verbatim-code>

  The first line must be present if we want to include quark mass effects. If
  the second line is present, top mass effects are included. If the third
  line is present, also bottom mass effects are included. If the fourth line
  is present (and also the third line is present), bottom mass effects are
  assumed to affect also the MiNLO Sudakov form factor, as discussed in
  <cite|masseffects>.

  In the setup of the <with|font-family|tt|PaperRun> directory, a first
  sample is generated with no quark mass effects. After that, the mass
  effects when only top is included, when both top and bottom are included,
  and when both top and bottom are included, including the effects of the
  bottom mass in the Sudakov form factor, are all added using the POWHEG
  reweighting facility.

  The <math|m<rsub|t>=\<infty\>> sample is generated with the following Les
  Houches labels:

  <\verbatim-code>
    lhrwgt_group_name 'MiNLO'

    lhrwgt_id 'mtinf'

    lhrwgt_descr 'mt infinity result'
  </verbatim-code>

  and no quark-mass flags set. This results in Les Houches event files
  including in the header the follwing block:

  <\verbatim-code>
    \<less\>header\<gtr\>

    \<less\>initrwgt\<gtr\>

    \<less\>weightgroup name='MiNLO'\<gtr\>

    \<less\>weight id='mtinf'\<gtr\> mt infinity result
    \<less\>/weight\<gtr\>

    \<less\>/weightgroup\<gtr\>

    \<less\>/initrwgt\<gtr\>

    \<less\>/header\<gtr\>

    \;
  </verbatim-code>

  \ The program is then run in reweighting mode, with the following variables
  settings:

  <\verbatim-code>
    compute_rwgt 1

    quarkmasseffects 1 \ \ ! Include quark mass effects

    topmass 172.5d0 \ \ \ \ \ ! top quark mass (172.5 in
    HNNLO-patches/mdata.f)

    lhrwgt_group_name 'MiNLO'

    lhrwgt_id 'mt'

    lhrwgt_descr 'Finite mt only'
  </verbatim-code>

  This run takes as input the name of an existing Les Houches event file,
  <verbatim|\<less\>file\<gtr\>.lhe>, and creates a new file, named
  \ <verbatim|\<less\>file\<gtr\>-rwgt.lhe>, including the new weight. This
  file is then renamed <verbatim|\<less\>file\<gtr\>.lhe>, and the procedure
  is repeated, with the following settings:

  <\verbatim-code>
    compute_rwgt 1

    quarkmasseffects 1 \ \ ! Include quark mass effects

    topmass 172.5d0 \ \ \ \ \ ! top quark mass (172.5 in
    HNNLO-patches/mdata.f)

    bottommass 3.38d0 \ \ \ ! bottom quark mass in MSbar at MH

    lhrwgt_group_name 'MiNLO'

    lhrwgt_id 'mtmb'

    lhrwgt_descr 'Finite mt and mb'
  </verbatim-code>

  and then again with the settings:

  <\verbatim-code>
    compute_rwgt 1

    quarkmasseffects 1 \ \ ! Include quark mass effects

    topmass 172.5d0 \ \ \ \ \ ! top quark mass (172.5 in
    HNNLO-patches/mdata.f)

    bottommass 3.38d0 \ \ \ ! bottom quark mass in MSbar at MH

    bmass_in_minlo 1

    lhrwgt_group_name 'MiNLO'

    lhrwgt_id 'mtmb'

    lhrwgt_descr 'Finite mt and mb, b mass in minlo'
  </verbatim-code>

  At the end of the procedure, the Les Houches file carries the following
  reweighting information in the header:

  <\verbatim-code>
    \<less\>header\<gtr\>

    \<less\>initrwgt\<gtr\>

    \<less\>weightgroup name='MiNLO'\<gtr\>

    \<less\>weight id='mtinf'\<gtr\> mt infinity result
    \<less\>/weight\<gtr\>

    \<less\>weight id='mt'\<gtr\> Finite mt only \<less\>/weight\<gtr\>

    \<less\>weight id='mtmb'\<gtr\> Finite mt and mb \<less\>/weight\<gtr\>

    \<less\>weight id='mtmb-bminlo'\<gtr\> Finite mt and mb, b mass in minlo
    \<less\>/weight\<gtr\>

    \<less\>/weightgroup\<gtr\>

    \<less\>/initrwgt\<gtr\>

    \<less\>/header\<gtr\>
  </verbatim-code>

  and after each event, we find reweighting information of the form:

  <\verbatim-code>
    \<less\>rwgt\<gtr\>

    \<less\>wgt id='mtinf'\<gtr\> \ 0.141990268E+02 \<less\>/wgt\<gtr\>

    \<less\>wgt id='mt'\<gtr\> \ 0.838351860E+01 \<less\>/wgt\<gtr\>

    \<less\>wgt id='mtmb'\<gtr\> \ 0.846010435E+01 \<less\>/wgt\<gtr\>

    \<less\>wgt id='mtmb-bminlo'\<gtr\> \ 0.845939328E+01 \<less\>/wgt\<gtr\>

    \<less\>/rwgt\<gtr\>
  </verbatim-code>

  <section|Generation of the <verbatim|hnnlo> output for reweighting>

  The fundamental NNLO ingredient required to obtain NNLO accurate event
  samples is, for gluon fusion Higgs production, the Higgs boson's rapidity
  spectrum.<footnote|In general, for NNLOPS simulations one requires an NNLO
  accurate distribution for some Born kinematics. In this case the Born
  kinematics is fully specified in terms of the rapidity of the boson.> This
  is obtained using the <verbatim|hnnlo> program of <cite|Grazzini:2008tf>
  <cite|Catani:2007vq> <cite|Grazzini:2013mca>. We here give instructions on
  how to obtain such distributions, suitable for combination with the
  <verbatim|HJ-MiNLO> events via the <with|font-family|tt|nllopsreweighter>
  program.

  <\enumerate-numeric>
    <item>Make sure that the <verbatim|lhapdf >package is
    installed:<next-line><hlink|https://lhapdf.hepforge.org|>

    <item>Download <verbatim|hnnlo> from the url:<next-line><verbatim|$ wget>
    <hlink|http://theory.fi.infn.it/grazzini/codes/hnnlo-v2.0.tgz|>

    <item>Unpack the tarball in a convenient location:

    <\verbatim-code>
      $ tar -xzvf hnnlo-v2.0.tgz

      $ ls hnnlo-v2.0
    </verbatim-code>

    Under the parent directory <with|font-family|tt|hnnlo-v2.0> one should
    find subdirectories <with|font-family|tt|bin>, <with|font-family|tt|doc>,
    <with|font-family|tt|obj>, <with|font-family|tt|src>, and a
    <with|font-family|tt|makefile>.

    <item>Enter the <verbatim|hnnlo> parent directory<next-line><space|5mm>
    <with|font-family|tt|$ cd hnnlo-v2.0>

    <item>Replace the <verbatim|hnnlo> default <with|font-family|tt|makefile>
    with the one from the <with|font-family|tt|HJ/NNLOPS-mass-effects>
    directory<next-line><space|5mm> <with|font-family|tt|$ cp
    /path/to/HJ/NNLOPS-mass-effects/HNNLO-makefile ./makefile>

    <item>Copy the <verbatim|hnnlo> patches directory into the parent
    directory<next-line><space|5mm> <with|font-family|tt|$ cp -r
    /path/to/HJ/NNLOPS-mass-effects/HNNLO-patches ./>

    <item>Build the code<next-line><space|5mm><with|font-family|tt|$
    make><next-line>A message <with|font-family|tt|<emdash>-\<gtr\> HNNLO
    compiled with LHAPDF routines \<less\><emdash>-> indicates success.

    <item>Copy over the required input file to the <verbatim|bin> directory.
    In the <verbatim|PaperRun> example:

    <space|5mm><with|font-family|tt|$ cd bin><next-line><space|5mm><with|font-family|tt|$
    cp ...POWHEG-BOX-V2/HJ/PaperRun/HNNLO*.input .>

    <item>Perform the run. For example:

    <\verbatim-code>
      $ ./hnnlo \<less\> HNNLO-LHC8-R04-APX0-11.input \<gtr\><space|0.17em>
      my.log
    </verbatim-code>
  </enumerate-numeric>

  This generates the <verbatim|HNNLO-LHC8-R04-APX0-11-10103.top> file. The
  file name reminds us that the run was performed with central
  renormalization and factorization scales (the 11 string),
  <verbatim|approxim 0>, and <verbatim|rseed 01103>. In order to obtain the
  <verbatim|HNNLO-LHC8-R04-APX0-11.top> file, 350 runs were performed. Each
  run takes of the order of 20 hours. The large number of runs is needed if
  one wants to get precise results also near the end-points of the rapidity
  spectrum.

  The output of the different runs were combined using the
  <verbatim|mergedata.f> program (located in the <verbatim|HJ/NNLOPS>
  directory). It is used as follows:

  <\verbatim-code>
    $ gfortran -o mergedata nergedata.f

    $ ./mergedata

    1

    filename1

    filename2

    ...
  </verbatim-code>

  where the list of filenames is terminated by an empty line.

  <section|Generation of NNLOPS reweighted events>

  In order to generate an NNLOPS sample, one needs first to compute certain
  distributions using the hnnlo program of <cite|Grazzini:2008tf>
  <cite|Catani:2007vq> <cite|Grazzini:2013mca>. This is documented in the
  HJ/NNLOPS/hnnlops.pdf document. In the <verbatim|PaperRun> directory, there
  are three input files and three output files of the <verbatim|hnnlo>
  program, with matching names. They correspond to the <verbatim|hnnlo> runs
  for the <math|m<rsub|t>=\<infty\>>, <math|m<rsub|t>> only, and
  <math|m<rsub|t>,m<rsub|b>> effects included. The data files (ending in
  <verbatim|.top>) are needed, together with the <verbatim|.lhe> file
  generated as described in the previous section, for the generation of the
  NNLOPS weights. These are obtained with the <verbatim|nnlopsreweighter>
  program, that is generated in the <verbatim|HJ> directory by simply doing

  <\verbatim-code>
    $ make nnlopsreweighter
  </verbatim-code>

  The <verbatim|nnlopsreweighter> program reads the input file
  <verbatim|nnlopsreweighter.input>, that contains the following lines:

  <\verbatim-code>
    lhfile pwgevents.lhe
  </verbatim-code>

  that instructs the program to read the events from the
  <verbatim|pwgevents.lhe> file. Then the section

  <\verbatim-code>
    nnlofiles

    'nn-mtinf' HNNLO-LHC8-R04-APX0-11.top

    'nn-mt' HNNLO-LHC8-R04-APX1-11.top

    'nn-mtmb' HNNLO-LHC8-R04-APX2-11.top
  </verbatim-code>

  specify the <verbatim|hnnlo> output file that should be read. Finally, we
  have the section

  <\verbatim-code>
    \<less\>initrwgt\<gtr\>

    \<less\>weightgroup name='nnl'\<gtr\>\ 

    \<less\>weight id='nnlops-mt'\<gtr\> combine 'nn-mtinf' and 'mtinf'
    \<less\>/weight\<gtr\>\ 

    \<less\>weight id='nnlops-mt'\<gtr\> combine 'nn-mt' and 'mt'
    \<less\>/weight\<gtr\>\ 

    \<less\>weight id='nnlops-mtmb'\<gtr\> combine 'nn-mtmb' and 'mtmb'
    \<less\>/weight\<gtr\>\ 

    \<less\>weight id='nnlops-mtmb-bminlo'\<gtr\> combine 'nn-mtmb' and
    'mtmb-bminlo' \<less\>/weight\<gtr\>\ 

    \<less\>/weightgroup\<gtr\>

    \<less\>/initrwgt\<gtr\>
  </verbatim-code>

  This has the format of a Les Houches header for reweighting. By running the
  <verbatim|nnlopsreweighter> program, we produce a new event file, named
  <verbatim|pwgevents.lhe.nnlo>, containing the following header

  <\verbatim-code>
    \<less\>header\<gtr\>

    \<less\>initrwgt\<gtr\>

    \<less\>weightgroup name='MiNLO'\<gtr\>

    \<less\>weight id='mtinf'\<gtr\> mt infinity result
    \<less\>/weight\<gtr\>

    \<less\>weight id='mt'\<gtr\> Finite mt only \<less\>/weight\<gtr\>

    \<less\>weight id='mtmb'\<gtr\> Finite mt and mb \<less\>/weight\<gtr\>

    \<less\>weight id='mtmb-bminlo'\<gtr\> Finite mt and mb, b mass in minlo
    \<less\>/weight\<gtr\>

    \<less\>/weightgroup\<gtr\>

    \;

    \<less\>weightgroup name='nnl'\<gtr\>

    \<less\>weight id='nnlops-mt'\<gtr\> combines 'nn-mtinf' and 'mtinf'
    \<less\>/weight\<gtr\>

    \<less\>weight id='nnlops-mt'\<gtr\> combines 'nn-mt' and 'mt'
    \<less\>/weight\<gtr\>

    \<less\>weight id='nnlops-mtmb'\<gtr\> combines 'nn-mtmb' and 'mtmb'
    \<less\>/weight\<gtr\>

    \<less\>weight id='nnlops-mtmb-bminlo'\<gtr\> combines 'nn-mtmb' and
    'mtmb-bminlo' \<less\>/weight\<gtr\>

    \<less\>/weightgroup\<gtr\>

    \;

    \<less\>/initrwgt\<gtr\>

    \<less\>/header\<gtr\>
  </verbatim-code>

  where the new entry mentioned in the <verbatim|nnlopsreweighter.input>
  appear, as well as the corresponding weights at the end of each event. The
  weights are obtained as follows. The <verbatim|nnlopsreweighter> program
  looks in the description line of each new weight mentioned in the
  <verbatim|nnlopsreweighter.input> file, looking for a string that appears
  as a weight id in the <verbatim|pwgevents.lhe> file, and for a string that
  appears as a prefix to a file name in the <verbatim|nnlofiles> section of
  the <verbatim|nnlopsreweighter.input> file. The new weight is then defined
  as the combination of the corresponding weight in the event file, and in
  the corresponding <verbatim|hnnlo> output file. In this way, any
  combination can be produced, and it is up to the user to use a sensible
  one. In the example reported in the <verbatim|PaperRun> directory, we
  combine <verbatim|HJ-MiNLO> and <verbatim|hnnlo> results that match in the
  accuracy of the quark mass effects that have been included. Note that the
  <verbatim|'mtmb'> and the <verbatim|'mtmb-bminlo'> weights are reweighted
  with the same <verbatim|'nn-mtmb'> tagged file, i.e. the
  <verbatim|HNNLO-LHC8-R04-APX2-11.top> file. In fact, there is no analog of
  the<verbatim| bmass_in_minlo> option in the <verbatim|hnnlo> program, that
  we only use to compute <math|p<rsub|T>> integrated distributions, without
  thus including any Sudakov form factor.

  <\bibliography|bib|JHEP.bst|paper.bib>
    <with|href|<macro|1|2|<arg|2>>|<\bib-list|1>
      <bibitem*|1><label|bib-masseffects>K.<nbsp>Hamilton, P.<nbsp>Nason, and
      G.<nbsp>Zanderighi, <with|font-shape|italic|Finite quark-mass effects
      in the NNLOPS POWHEG+MiNLO Higgs generator>,
      <hlink|<with|font-family|tt|xxxx.xxxx>|http://xxx.lanl.gov/abs/xxxx.xxxx>.

      <bibitem*|2><label|bib-Hamilton:2013fea>K.<nbsp>Hamilton,
      P.<nbsp>Nason, E.<nbsp>Re, and G.<nbsp>Zanderighi,
      <with|font-shape|italic|NNLOPS simulation of Higgs boson production>,
      <with|font-shape|italic|JHEP> <with|font-series|bold|1310> (2013) 222,
      [<hlink|<with|font-family|tt|1309.0017>|http://xxx.lanl.gov/abs/1309.0017>].

      <bibitem*|3><label|bib-Hamilton:2012rf>K.<nbsp>Hamilton, P.<nbsp>Nason,
      C.<nbsp>Oleari, and G.<nbsp>Zanderighi, <with|font-shape|italic|Merging
      H/W/Z + 0 and 1 jet at NLO with no merging scale: a path to parton
      shower + NNLO matching>, <hlink|<with|font-family|tt|1212.4504>|http://xxx.lanl.gov/abs/1212.4504>.

      <bibitem*|4><label|bib-Grazzini:2008tf>M.<nbsp>Grazzini,
      <with|font-shape|italic|NNLO predictions for the Higgs boson signal in
      the H <math|\<rightarrow\>> <math|W<rsup|+>*W<rsup|->>
      <math|\<rightarrow\>> <math|<wide|l|\<bar\>>\<nu\>>
      <math|l<wide|\<nu\>|\<bar\>>> and H <math|\<rightarrow\>> ZZ
      <math|\<rightarrow\>> 4l decay channels>, <with|font-shape|italic|JHEP>
      <with|font-series|bold|0802> (2008) 043,
      [<hlink|<with|font-family|tt|0801.3232>|http://xxx.lanl.gov/abs/0801.3232>].

      <bibitem*|5><label|bib-Catani:2007vq>S.<nbsp>Catani and
      M.<nbsp>Grazzini, <with|font-shape|italic|An NNLO subtraction formalism
      in hadron collisions and its application to Higgs boson production at
      the LHC>, <with|font-shape|italic|Phys.Rev.Lett.>
      <with|font-series|bold|98> (2007) 222002,
      [<hlink|<with|font-family|tt|hep-ph/0703012>|http://xxx.lanl.gov/abs/hep-ph/0703012>].

      <bibitem*|6><label|bib-Grazzini:2013mca>M.<nbsp>Grazzini and
      H.<nbsp>Sargsyan, <with|font-shape|italic|Heavy-quark mass effects in
      Higgs boson production at the LHC>, <with|font-shape|italic|JHEP>
      <with|font-series|bold|1309> (2013) 129,
      [<hlink|<with|font-family|tt|1306.4581>|http://xxx.lanl.gov/abs/1306.4581>].
    </bib-list>>
  </bibliography>
</body>

<initial|<\collection>
</collection>>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-2|<tuple|2|?>>
    <associate|auto-3|<tuple|3|?>>
    <associate|auto-4|<tuple|3|?>>
    <associate|bib-Campbell:2012am|<tuple|1|?>>
    <associate|bib-Campbell:2013vha|<tuple|2|?>>
    <associate|bib-Catani:2007vq|<tuple|5|?>>
    <associate|bib-Grazzini:2008tf|<tuple|4|?>>
    <associate|bib-Grazzini:2013mca|<tuple|6|?>>
    <associate|bib-Hamilton:2012rf|<tuple|3|?>>
    <associate|bib-Hamilton:2013fea|<tuple|2|?>>
    <associate|bib-masseffects|<tuple|1|?>>
    <associate|footnote-1|<tuple|1|?>>
    <associate|footnote-2|<tuple|2|?>>
    <associate|footnote-3|<tuple|3|?>>
    <associate|footnr-1|<tuple|1|?>>
    <associate|footnr-3|<tuple|3|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      masseffects

      Hamilton:2013fea

      Hamilton:2012rf

      masseffects

      masseffects

      Grazzini:2008tf

      Catani:2007vq

      Grazzini:2013mca

      Grazzini:2008tf

      Catani:2007vq

      Grazzini:2013mca
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Production
      of the HJ-MiNLO events with mass effects.>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Generation
      of the <with|font-family|<quote|tt>|language|<quote|verbatim>|hnnlo>
      output for reweighting> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Generation
      of NNLOPS reweighted events> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>