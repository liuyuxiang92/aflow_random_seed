****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2020                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
****************************************************************************************************
LATEST VERSION OF THE FILE:                                   materials.duke.edu/AFLOW/README_PROTO.TXT
****************************************************************************************************
****************************************************************************************************
****************************************************************************************************
PROTOTYPES LIB1

STR#           Symmetry  Name          TYPE      PROTO       SG       conc
-----------------------------------------------------------------------------------------------
FCC
A1,FCC,1,2        fcc/cF4   Fm-3m         A1        Cu          #225
A2,BCC,58,59      bcc/cI2   Im-3m         A2        W           #229
A3,HCP,115,117    hcp/hcp2  P6_3/mmc      A3        Mg          #194
A4,DIA,301,302    cF8       Fd-3m         A4        C (DIA)     #227  books
A5,betaSn,305,306 tI4       I41/amd       A5        betaSn      #141  Navy Database
A6,In,303,304     tI2       I4/mmm        A6        Indium      #139  Navy Database
A7,As,307,308     hR6       R-3m          A7        alphaAs     #166  Kolmogorov
A8
A9
A10,Hg,349,350    hR1       R-3m          A10       alphaHg     #166  ICSD

****************************************************************************************************
****************************************************************************************************
****************************************************************************************************
PROTOTYPES BINARY

STR#           Pearson   Symmetry      Name     TYPE   PROTO      SG      notes
-----------------------------------------------------------------------------------------------
FCC
1,2            fcc/cF4   Fm-3m         A1        Cu              #225                              G(pic1N1)
1.cub,2.cub    fcc/cF4   Fm-3m         A1        Cu              #225                              G(pic1N1)
3              tP2       P4/mmm        L1_0      AuCu            #123     AB  along [001]          G(pic1N2,1)
4              hR2       R-3m          L1_1      CuPt            #166     AB  along [111]          G(pic1N2,2)
5,6            tI6       I4/mmm        beta1/2   AB2             #139     AB2 along [100]z         G(pic1N3,2,"C11_b")
7,8            oI6       Immm          -         MoPt2           #71      AB2 along [011]          G(pic1N3,1)
9,10           hP3       P-3m1         alpha1/2  AB2             #164     AB2 along [111]z         G(pic1N3,3,"C6")
11,12          oS8       Cmmm          "L1_3"(unofficial) AB3    #65      "(A.5B.5)B along [111]"  M(CdPt3,PdPt3,Pd3Pt) G(pic1N4,2)
13,15          tP4       P4/mmm        Z1/3      AB3             #123     AB3 along [001]z         M(PtTc3)  G(pic1N4,8)
14             tP4       P4/nmm        Z2        A2B2            #129     A2B2 along [001]z        G(pic1N4,3,"B11")
16,18          mS8       C2/m          ?         AB3             #12      AB3 along [113]          G(pic1N4,12)
17             mS8       C2/m          W2        A2B2            #12      A2B2 along [311]z        G(pic1N4,7)
19,21          oP4       Pmmm          Y1/3      AB3             #47      AB3  along [011]z        G(pic1N4,9)
20             oP4       Pmmm          Y2        A2B2            #47      A2B2 along [011]z        M(AgZr) G(pic1N4,4)
22,24          tI8       I4/mmm        D0_22     Al3Ti           #139     AB3  along [201]         G(pic1N4,11)
23             tI8       I41/amd       CH "40"   NbP             #141     A2B2 along [201]         G(pic1N4,6,NbP)
25,26          cP4       Pm3(bar)m     L1_2      AuCu3           #221     "(A.5B.5)B  along [001]" G(pic1N4,1)
27,29          hR4       R-3m          V1/3      AB3             #166     AB3  along [111]z        G(pic1N4,10)
28             hR4       R-3m          V2        A2B2            #166     A2B2 along [111]z        G(pic1N4,5)
BCC
58,59          bcc/cI2   Im-3m         A2        W               #229
58.cub,59.cub  bcc/cI2   Im-3m         A2        W               #229
60             oS4       Cmmm          D^19_2h   gamma-IrV       #65      AB along [101]
61             cP2       Pm-3m         B2        CsCl            #221     AB along [001]
62,63          hP3       P-3m1         ?         AB2             #164     AB2 along [211]
64,65          oF12      Fmmm          ?         AB2             #69      AB2 along [011] WARNI FIX NpPd
66,67          tI6       I4/mmm        C11_b     MoSi2           #139     AB2 along [001] WARNI FIX NpPd
68,69          hR4       R-3m          ?         AB3             #166     UNDEF
70,72          oS8       Cmmm          ?         AB3             #65      UNDEF
71             oS8       Cmma          ?         A2B2            #67      UNDEF
73,75          mP4       P2/m          ?         AB3             #10      UNDEF
74             mP4       P21/m         ?         A2B2            #11      UNDEF
76,78          tP4       P4/mmm        ?         AB3             #123     AB3 along [001]z
77             tP4       P4/nmm        B11       gamma-CuTi      #129     A2B2 along [001]z
79,81          oI8       Immm          ?         AB3             #71      UNDEF  M(MoTi3, Mo3Ti, Nb3Tc,MoTi3, Mo3Ti, Nb3Tc)
80             oI8       Imma          ?         A2B2            #74      UNDEF  M(MoTi)
82,83          tP4       P4/mmm        L6_0      CuTi3           #123
84,86          cF16      Fm-3m         D0_3      AlFe3/BiF3      #225
85             cF16      Fd-3m         B32       NaTl            #227
HCP
115,117        hcp/hcp2  P6_3/mmc      A3        Mg              #194
116            hP2       P_6m2         B_h       WC              #187
118,121        oP4       Pmm2          ?         AB3             #25      UNDEF
119            oP4       Pmmn          -         CuTe            #59
120            oP4       Pmma          B19       AuCd            #51
123            oI8       C2/m          ?         A2B2            #12      UNDEF
122,124        oI8       Imm2          ?         AB3             #44      UNDEF
125,127        hP4       P-6m2         ?         AB3             #187     UNDEF
126            hP4       P-3m1         ?         A2B2            #164     UNDEF M(RhRu)
128,132        mS12      Cm            ?         AB5             #8       UNDEF
129,134        mS12      C2/m          ?         A4B2            #12      BiHf2-134       #12  Hf-Paper
130,137        mS12      C2/m          ?         A4B2            #12
131,135        mS12      C2/m          ?         A4B2            #12
133,140        mS12      Cm            ?         A2B4            #8       UNDEF
136            mS12      Cm            ?         A3B3            #8       UNDEF
138,139        mS12      Cm            ?         A3B3            #8       UNDEF
141,145        oS12      Amm2          ?         AB5             #38
142,147        oS12      Cmcm          ?         A4B2            #63      M(RhRu2)
143,150        oS12      Cmcm          ?         A4B2            #63      mistaken for C49, beware
*144,148       oS12      Cmcm          ?         A4B2            #63      NOT RUN
*146,153       oS12      Amm2          ?         A2B4            #38      NOT RUN
*149,151       oS12      Amm2          ?         A3B3            #38      NOT RUN
*152,157       oS12      Amm2          ?         A3B3            #38      NOT RUN
*154,160       hP6       P-62m         ?         AB5             #189     NOT RUN
*155,158       hP6       P-62m         ?         A2B4            #189     NOT RUN
*156,159       hP6       Cmcm          ?         A2B4            #63      NOT RUN
*161,167       hR6       R32           ?         AB5             #155     NOT RUN
*162,165       hR6       R32           ?         A2B4            #155     NOT RUN
*163,166       hR6       C2/c          ?         A2B4            #15      NOT RUN
*164           hR6       C2            ?         A3B3            #5       NOT RUN
*168,169       hP6       P-6m2         ?         AB5             #187     NOT RUN
*170,175       hP6       P-6m2         ?         A2B4            #187     NOT RUN
*172           mS8       P-6m2         ?         A3B3            #187     NOT RUN
*171,176       hP6       P-3m1         ?         A2B4            #164     NOT RUN
*173,177       hP6       P63/mmc       ?         A2B4            #194     NOT RUN
*174           hP6       P3m1          ?         A3B3            #156     NOT RUN
LAVES AND OMEGAS
178,179        hP12      P6_3/mm       C14       MgZn2           #194    A4B8  Hex Laves Str
180,181        oP16      Pnma          D0_11     Fe3C            #62
182,183        cF24      Fd-3m         C15       Cu2Mg           #227    Cub Laves Phase
184,185        cP8       Pm3(bar)n     A15       Cr3Si           #223
186,187        hP8       P6_3/mmc      D0_19     Ni3Sn           #194
188,189        oS12      Cmcm          C49       ZrSi2           #63
190,191        hP3       P3m1      omega phase with Z=1/4 !      #164
192,193        oS8       Cmcm          B33B_f    CrB             #63
194            oS8       Cmcm                    ITl ITl         #63     Pauling
195            cP8       P2_13         B20       FeSi            #198    Pauling
196,197        tP4       P4/nmm O2               CdTi            #129    Pauling
198            oP24      Pbcm                    InTh            #57     Pauling
199            tI16      I4/mcm        B37       SeTl            #140    Pauling  

201            cF8       Fm-3m         B1        NaCl            #225
202,203        tI10      I4/mmm        D1_3      BaAl4           #139    NEW
204,205        hP6       P6/mmm        D2_d      CaCu5           #191
206,207        tI26      I4/mmm        D2_b      ThMn12          #139
208,209        hP9       P-62m         C22       Fe2P            #189
210,211        oP12      Pnma          C37       Co2Si           #62     AL1851.cif
212,213        hR19                              Th2Zn17         #166
214,215        cI28      I-43d         D7_3      Th3P4           #220
216,217        hP3       P6/mmm        C32       AlB2            #191
218            cF8       F-43m         B3        ZnS             #216
219            hP4       P63/mc        B4        ZnS             #186
220,221        hP4       P63/mmc       B8_1      NiAs            #194
222,223        hP16      P63/mcm       D8_8      Mn5Si3          #193    128601.cif
224,225        hP38      P63/mmc                 Th2Ni17         #194
226,227        cP8       Pm-3n         A15       Cr3Si           #223
228,229        tP6       P4/nmm        C38       Cu2Sb           #129
230,231        tI12                    C16       CuAl2           #140
232,233        tI6       I4/mmm        C11_b     MoSi2           #139
234,235        tl12      I4/mcm        C16       Al2Cu           #140    130617.cif
236,237        hP142     P63/mmc                 Cd58Y13         #194   ICSD_102087
****NO 236,237        cl184     Im-3          Cd6Y     (Y12Cd68)        #204    130045.cif
****NO 238,239        hp24      P-6m2                   Co3V            #187    130340.cif
238,239        hp24      P63/mmc       Pu3Al                     #194    Pauling FILE = Co3V
***240,241               P42/mnm       D8_b      CrFe            #136    130307.cif ???? HARD
242,243        oP24      PmnmO2                  NbPd3           #59     130379.cif
244,245                  P63/mmc       D0_24     Ni3Ti           #194    AL5957.cif
246,247        hP6       P63/mmc                 CaIn2           #194    AL5296.cif
248,249        tI12      I41/amd       C_c       ThSi2           #141    AL2272.cif
***250,251     tI6       I4/mmm        C11_a     CaC2            #139    AL1909.cif UNRUNNABLE
252,253        cF24      F-43m         C15_b     AuBe5           #216    AL3669.cif
254,255        tI32      I4/mcm                  W5Si3           #140    AL2650.cif
               hP5       P63/mmc       D5_2      La2O3           #194    130882.cif FRAC OCCUP
256            oP8       Pnma          B27       BFe             #62     AL1744.cif
257,258        oP6       Pnnm          C18       FeS2            #58     marcasite AL1657.cif
259,260        hR5       R-3m          C33       Bi2Te3          #166    120146.cif
261,262        cF96      Fd-3m                   Ti2Ni           #227
263,264        tI14      I4/mmm                  Ti3Cu4          #139    120654.cif
265,266        tP10                              Ti2Cu3          #129    IMPOSSIBLE
267,268        oP6       Pnnm          C18       FeS2            #52     marcasite AL1657.cif
269,270        hP3       P-3m1         C6        CdI2            #164    AL4899.cif
271,272        oS16      Cmcm                    Cd3Y            #63     18260.cif
273,274        hP6       P63/mmc       B8_2      Ni2In           #194    28249.cif
273p,274p      hP5       P-6m2      Fe1.5Ge1 partial of Ni2In    #187    Pauling S525640
275,276        hP6       P63/m                   Ni2Si           #176    AL5234.cif
277,278        oP8       PmmnO1        D0_a      betaCu3Ti       #59     AL1729.cif
279,280        tI16      I4/mmm        D0_23     Al3Zr           #139    Navy Database
281,282        cP12      Pa-3          C2        FeS2            #205    chris pyrite
283,284        oI12      Imma                    GdSi2(1.4)      #74     chris z=(0.0362,0.2138,0.6250)
285,286        tI10      I4/m          D1_a      MoNi4           #87     Navy Database
287,288        tI6       I4/mmm                  CuZr2           #139    Pauling file
289,290        cP4       Pm-3m         D0_9      alphaReO3       #221    Navy Database
291,292        tP4       P4/nmm        B10       PbO/PbS         #129    Navy Database
301,302        cF8       Fd-3m         A4        C (DIA)         #227    books
303,304        tI2       I4/mmm        A6        Indium          #139    Navy Database
305,306        tI4       I41/amd       A5        betaSn          #141    Navy Database
307,308        hR6       R-3m          A7        alphaAs         #166    Kolmo
309,310        cF32      Fm-3m                   CuPt7           #225    (Ca7Ge) Navy Database
311,312        tI18      I4/mmm                  NbNi8           #139    (Pt8Ti) Navy Database
313,314        tI18      I4/mmm                  V4Zn5           #139    Navy Database
315,316        hP24      P6_3/mmc      C36       MgNi2           #194    Navy Database
317,318        oS8       Cmca          A11       alphaGa         #64     Navy Database
319,320        cI58      I(-4)3m       A12       alphaMn         #217    Navy Database
321,322        cP20      P4_132        A13       betaMn          #213    Navy Database
323,324        hP3       P3_21-D_3^6   A8        gammaSe         #152    Navy Database
325,326        hP4       P63/mmc       A9        graphite        #194    Navy Database
327,328        tI2       I4/mmm        Aa        alphaPa         #139    Navy Database
329,330        tP30      P4_2/mnm      Ab        betaU           #136    Navy Database/Pauling
331,332        oP8       I4/mmm        Ac        alphaNp         #62     Navy Database/Pauling
333,334        tP4       P4/nmm        Ad        betaNp          #129    Navy Database/Pauling
335,336        hP1       P6/mmm        Af hex    HgSn_{6-10}     #191    Navy Database
337,338        todo Ag T-50 boron
339,340        todo Ah simple cubic
341,342        todo Ai betaPo
343,344        todo Ak alphaSe
345,346        todo Al betaSe
347,348        todo A3' alphaLa
349,350        hR1       R-3m          A10       alphaHg         #166    ICSD   A10.A/B
351,352        todo A14 In2
353,354        todo A16 alphaS
355,356        todo A17 black P
357,358        todo A20 alphaU
359,360        hP18     P63/mcm                  Ga4Ti5          #193    Pauling
361,362        tP20     P42/mnm                  Al2Zr3          #136    Pauling
363,364        hP7      P6/mmm                   Al3Zr4          #191    Pauling
365,366        oF40     F2dd                     Al3Zr2          #43     Pauling
367,368        cF112    Fm-3c                    NaZn13          #226    Pauling
369,370        cP7      Pm-3m                    CaB6            #221    Pauling
371,372        hP5      P-3m1          D5_19     Al3Ni2          #164    Navy Database
373,374        cI40     I-43m          Ga4Ni(Ga3.62Ni0.97)       #217    Navy Database (vacancies I23,     #197)
375,376        oS16     Cmmm                     Ga3Pt5          #65     Navy Database
377,378        mS44     C12/m1                   Ga9Ni13         #12     Pauling
379,380        cI112    Ia(-3)d                  Ga4Ni3          #230    Navy Database
381,382        tI24     I4_1/amd                 Ga2Hf           #141    Navy Database
383,384        cF116    Fm-m                     C6Mn23          #225    ICSD
385,386        oP44     Pnma                     Cu8Hf3          #62     ICSD
387,388        oS68     Cmca                     Ni10Zr7         #64     ICSD_240191
389,390        tP3      P4/mmm                   Hg2Pt           #123    ICSD_104339
391,392        oP20     Pnma                     B4Mg            #62     Pauling
393,394        oI64     Imma                     B7Mg            #74     Pauling
395,396        oP16     Pbam                     Ge3Rh5          #55     Pauling
397,398        tP14     P4/mbm                   In4Ti3          #127    Pauling
399,400        tP34     P4/nmm.O2      V8Sb9     V7.49Sb9        #129    Pauling
401,402        tI160    I-42d                    Li3B14          #122    Pauling/Kolmogorov
403,404        tP20     Pa/mbm                   LiB3            #127    Pauling/Kolmogorov
405,406        hR4      R-3m                     LiB-MS1         #166    Pauling  
407,408        hP8      P63/mmc                  LiB-MS2         #194    Pauling    

411,412        cI58     I-43m                    Re24Ti5         #217    Pauling
413,414        mS36     C12/m1                   Ni7Zr2          #12     Pauling
415,416        tI40     I4/m                     Pt11Zr9         #87     Pauling
417,418        oP20     Pnma                     Au4Zr           #62     Pauling
419,420        oI142    Immm                     Hf54Os17        #71     Pauling
421,422        hP9      P6_222         C40       CrSi2           #180    Pauling
423.424.t0     hP68     P6/m                  Ag51Gd14.o3.t0     #175    Pauling
423.424.t1     hP68     P6/m                  Ag51Gd14.o3.t1     #175    Pauling
423.424.t2     hP68     P6/m                  Ag51Gd14.o3.t2     #175    Pauling
423,424.t3     hP68     P6/m                  Ag51Gd14.o2.t3     #175    Pauling
423,424.t4     hP68     P6/m                  Ag51Gd14.o2.t4     #175    Pauling
423,424.t5     hP68     P6/m                  Ag51Gd14.o2.t5     #175    Pauling
423,424.t6     hP68     P6/m                  Ag51Gd14.o2.t6     #175    Pauling
423,424.t7     hP68     P6/m                  Ag51Gd14.o1.t7     #175    Pauling
423,424.t8     hP68     P6/m                  Ag51Gd14.o1.t8     #175    Pauling
423,424.t9     hP68     P6/m                  Ag51Gd14.o1.t9     #175    Pauling
425,426        cF184    Fd-3m                    Zn22Zr          #227    Pauling
427,428        oP24     Pnnm                     As2Ti           #58     Pauling
429,430        tI84     I4/mmm                   Ge10Ho11        #139    Pauling
431,432        hP48     P6_122                   Ir3Zr5          #178    Pauling
433,434        oP6      Pnnm                     CaCl2           #58     Pauling
435,436        cP5      P-43m                    CFe4            #215    Pauling
437,438        mS28     C12/c1                   C2Mn5           #15     Pauling
439,440        oP40     Pnma                     C3Mn7           #62     Pauling
441,442        hP20     P6_3mc                   Fe3Th7          #186    Pauling
443,444        hR24     R-3cr                    F3Fe            #167    Pauling
445,446        hR24     P321                     F3Fe            #150    Pauling  
447,448        mP4      P12_1/m1                 NiTi            #11     Pauling
449,450        cP6      Pn-3m          C3        Ag2O            #224    ICSD_281041 Cuprite
451,452        hP9      P-3            Bb        etaAg2Zn        #147    Navy Database
453,454        mS16     C2/m           D0_15     AlCl3           #12     ICSD_39566
455,456        hR8      R32            D0_14     AlF3            #155    ICSD_30274
457,458        hR10     R-3c           D5_1      alphaAl2O3      #167    ICSD_75479 Corundum
459,460        cI26     Im-3                     Al12W           #204    ICSD_58207
461,462        oP12     P212121                  Ag2Se           #19     ICSD_15213 Naumannite
463,464        hR26     R-3                      AlPd            #148    ICSD_58112
465,466        cP16     Pa-3           SC16      AlSb            #205    ICSD_41994
467,468        cI52     I-43m         gammaBrass Cu5Zn8          #217    Pauling
469,470        hR36     R-3m h                   BaPb3           #166    Pauling
471,472        hP8      Cmcm                     Hf3Sc*-h321     #63     Hf-Paper
473,474        hP6      P-62m                    Hf5Sc-h51       #189    Hf-Paper
475,476        mS12     C2/m                     BiHf2-134       #12     Hf-Paper
477,478        tP6      P4/mmm                   Hf5Pb-f63       #123    Hf-Paper
479,480        oS12     Cmmm                     HfPd5           #65     Hf-Paper
481,482        hR276    R-3ch                    Re25Zr21        #167.2  Pauling
483,484        tP6      P4/mmm        B11_3                      #123    Mike Mehl bcc_AAABBB^[001]
485,486        tP6      P4mm          B11_3'                     #99     Mike Mehl bcc_AABBAB^[001]
487,488        tP6      P4/mmm        Z3                         #123    Mike Mehl fcc_AAABBB^[001]
489,490        tP6      P4mm          Z3'                        #99     Mike Mehl fcc_AABBAB^[001]
491,492        oP102    Pmn2_1        Al13Co4                    #31     Pauling - Grin J_Al_Comp_206_243_1994
493,494        ms102    C1m1          Al12.13Co4   eta           #8      Pauling - Hurd Acta Cryst 15 441 1962
495,496        hP28     P6_3/mmc      Al5Co2     D8_11           #194    Pauling
497,498        mp22     P12_1/c1      Al9Co2                     #14     Pauling
499,500        tP10     P-42_1c                  Pd4Se           #114    Pauling
501,502        cF1192   Fd(-3)m                  NaCd2           #227    Samson, Nature 195, 259 (1962)
503,504        hP10     P6_3/m                   Er3Ru2          #176    Pauling
505,506        oP40     Pnma                     Pt3Sr7          #62     Pauling
507,508        cF120    Fm-3m                    Ir4Sc11         #225    Pauling
509,510        oP276    Pnna                     Ru25Y44         #52     Pauling
511,512        cP252    P4_132                   RuZn6           #213    Pauling
513,514        mP4      P21/m         B19'       NbRu-Beta''     #11     Eur.Phys.J.B 72 575(2009)
515,516        cF96     Fd-3m                    delta(CdNi)     #227.2  Pauling
517,518        hP3      P-3m1                    Cd2Ce           #164    Pauling
519,520        hP3      P6/mmm                   Hg2U            #191    Pauling
521,522        hP9      P-62m                    InMg2           #189    Pauling
523,524        oS16     Cmcm                     Cd3Er           #63     Pauling
525,526        hP8      P6_3/mmc                 CdMg3           #194    Pauling
527,528        hP24     P6_3/mmc                 CeNi3           #194    Pauling
529,530        oP32     Pnma                     CoSc3           #62     Pauling
531,532        tP4      P4/mmm                   Pb3Sr           #123    Pauling
533,534        hR36     R-3mh                    PuNi3           #166    Pauling
535,536        oP16     Pnma                     YZn3            #62     Pauling
537,538        hR4      R-3m                     Co2Y2*          #166    hcp-paper definition
539,540        mS12     C2/m                     Sc2Zr*          #12     hcp-paper definition
541,542        oI8      Immm                     Mo3Ti*-81       #71     hcp-paper/monster
543,544        oI8      Imma                     MoTi*-80        #74     hcp-paper/monster
545,546        hR3      P-3m1                    ReTi2*-81       #164    hcp-paper
547,548        tI6      I4/mmm                   Hf2Tl*-6        #139    hf-paper
549,550        oF12     Fmmm                     Be2Zn*-65       #69     hcp-paper
551,552        oI8      Imm2                     Re3Ru*-124      #44     hcp-paper
553,554        hP24     P63cm          D0_21     Cu2.82P         #185    Pauling
555,556        oP108    Pnma                     Mg2Au           #62     Pauling
557,558        cF16     Fm-3m                    BiF3            #225    Pearson's    crystal data
559,560        cI58     I-43m                    Al12Mg17        #217    Pauling
561,562        oS160    Cmcm                     MgAu3-x         #63     Pauling
563,564        oS64     Cmcm                     MgAu3+x         #63     Pauling
565,566        oF48     Fddd                     Mg2Cu           #70     Pauling
567,568        hP18     P-62c                    Mg2Ga           #190    Pauling
569,570        oP24     Pbam                     MgGa2           #55     Pauling
571,572        oI28     Ibam                     Mg5Ga2          #72     Pauling
573,574        tI28     I4/mmm                   Mg2Ga5          #139    Pauling
575,576        tI32     I41/a                    MgGa            #88     Pauling
577,578        hR72/24  R32                      Mg3Hg           #155    Pauling
579,580        hR48/16  R-3m                     Mg3In           #166    Pauling
581,582        cF408    F-43m                    Mg44Rh7         #216    Pauling
583,584        cF12     Fm-3m                    CaF2            #225    Pauling
585,586        hP8      P63/mmc        D0_18     Al3Ir           #194    Pauling
587,588        cF116    Fm-3m                    Mn23Th6         #225    Pauling
589,590        cP39     Pm-3                     Mg2Zn11         #200    Pauling
591,592        mS110    C2/m                     Mg4Zn7          #12     Pauling
593,594        cP20     P4_132                   Mg3Ru2          #213    Pauling
595,596        oS12     Cmcm                     Au2V            #63     Pauling
597,598        hP94     P63/mmc                  Sr9Mg38         #194    Pauling
600.ABCDE      tP30     P4_2/mnm          sigma  CrFe            #136    Pauling
600.XXXXX where X can be A,B (a total of 32 poscars)
611,612        cP140    Pm-3                     Rh13Sc57        #200    Pauling
613,614        hR13     Imm2                     Fe7W6           #166    Pauling
615,616        oP24     Pbca                     AuSn2           #61     Pauling
617,618        mS12     C2./m                    Bi2Pd           #12     Pauling
619,620        mS28     C2./m                    Bi2Pd5          #12     Pauling
621,622        oP16     Pnma                     Bi3Ni           #62     Pauling
623,624        oP32     Pnma                     Bi3Y5           #62     Pauling
625,626        cI120    Ia-3d                    Bi4Rh           #230    Pauling
627,628        hR48     R-3m                     BiMn3           #166    Pauling
629,630        oP16     Pmma                     BiPd3           #51     Pauling
631,632        tI12     I4/mmm                   BiTi2           #139    Pauling
633,634        mP12     P2_{1}./c                CoSb2           #14     Pauling
635,636        hR78     R3m                      Cu7Hg6          #160    Pearson (Gamma-Brass Cu5Zn8 Modified)
637,638        oI12     Imma                     Hg2K            #74     Pauling
639,640        cI52     I-43m                    Ir2Zn11         #217    Pauling files (Gamma-Brass Cu5Zn8 Modified Occupancy
641,642        tI32     I-4            D0_e      Ni3P            #82     Pauling
643.AB/BA      aP32     P-1                      KP15            #2      Pauling
644.AB/BA      tI116    I4_1/a                   Pu28Zr          #88     Pauling
645.AB/BA      cF52     Fm-3m                    B12U            #225    Pauling
646.AB/BA      oS68     Cmcm                     TiZn16          #63     Pauling
647.AB/BA      tI140    I4/mcm                   Rh2Y3           #140    Pauling
648.AB/BA      tP32     P4/ncc                   Rh3Pu5          #130    Pauling
649.AB/BA      tI32     I-42m                    alpha-SV3       #121    Pauling
650.AB         oP12     Pmma                     IrTa            #51     Pauling
651.AB/BA      oP36     Pnma                     Ge4Sm5          #62     Pauling
652.AB/BA      oP12     Pnma                     Cl2Pb           #62     Pauling
653.AB/BA      mP16     P12_1                    NbPt3           #11     Pauling
654.AB/BA      cI10     Im-3m                    Hg4Pt           #229    Pauling
655.AB/BA      hR14     R-3                      Pd4Pu3          #148    Pauling
656.AB/BA      hR15     R-3                      Er3Ni2          #148    Pauling
657.AB/BA      oI6      Immm                     Pd2Ti           #71     Pauling
658.AB/BA      tP8      P4/mmm                   Pd5Ti3          #123    Pauling
659.AB/BA      oS20     Cmcm                     Pd3Ti2          #63     Pauling
660.AB/BA      tP14     P4/mbm                   Hg5Mn2          #127    Pauling
661.AB/BA      oS32     Cmcm                     Pu3Pd5          #63     Pauling
662.AB/BA      oS12     Cmmm                     Ga2Zr           #65     G. Ghosh and M.Asta, Acta Mat. 53, 3225-3252 (2005).
663.AB/BA      mP30     P21/c                    FeP4_mP30       #14     Pauling
664.AB/BA      mS40     C2/C                     FeP4_mS40       #15     Pauling
665.AB/BA      tP16     P4/2mnm                  In3Ir           #136    Pauling
666.AB/BA      mS44     C2/m                     Fe6Ge5          #12     Pauling
667.AB/BA      mS42     C2/m                     Fe3Ga4          #12     Pauling
668.AB/BA      oS28     Cmcm                     Al6Mn           #63     Pauling
669.AB         oP8      Pnma                     AsFe            #62     Pauling
670.AB/BA      hP6      P6/mmm      B35          CoSn            #191    Pauling
671.AB/BA      cI32     Im-3        D0_2         As3Co           #204    Pauling
672.AB/BA      cI32     Im-3m                                    #229    Mike - 16 atom scell fcc, A15B1
673.AB/BA      cP16     Pm-3m                                    #221    Mike - 16 atom scell bcc, A15B1
674.A/B        hP3      P6/mmm                                   #191    Mike - omega phase (C32) in pure form
675.A/B        cP8      Pm-3n                                    #223    Mike - A15 phase in pure form
676.AB/BA      hP22     P6/mmc        Fe6.5Ge4-ideal(Fe14Ge8)    #194    Pauling
676p.AB/BA     oS42     Amm2          Fe6.5Ge4-part(Fe13Ge8)     #38     Pauling
677.AB/BA      hR26     R3m                      Al8Cr5          #160.2  ICSD_606753
678.AB/BA      tI24     I4mm                     Co5Ge7          #107    ICSD_43386
679.AB/BA      cP5      Pm-3m         E2_1                       #221    NAVY
680.AB/BA      tP12     P4/mmm                   beta-NW2        #123    NAVY (Hagg, sort of)
681.AB/BA      hR3      R3m                      alpha-MoS2      #160    NAVY
682.AB/BA      hR13     R-3m          D8_5r                      #166    NAVY
683.AB/BA      cI50     Im-3m                    Na13Cl12v7      #229    NAVY
684.A/B        cP8      P2_13                    alpha-N2        #198    NAVY
685.A/B        cP8      Pa-3                     alpha-N2p       #205    NAVY
686.A/B        hP4      P6_3/mmc                 beta-N2         #194    NAVY
687.A/B        hR16     R-3c                     epsilon-N2      #167    NAVY
688.A/B        tP4      P4_2/mnm                 gamma-N2        #136    NAVY
689.AB/BA      hP6      P6_3/mmc                 WN2_194af       #194    NAVY
690.AB/BA      mP12     P2_1/c                   Ag2Te           #14     NAVY
691.AB/BA      cF12     Fm-3m         C1                         #225    NAVY
692.AB/BA      hR3      R-3m          C19                        #166    NAVY
693.AB/BA      tP6      P4_2/mnm      C4                         #136    NAVY
694.AB/BA      hP6      P6_3/mmc      C7                         #194    NAVY
695.AB/BA      hP3      P-6m2                    WN2_187fg       #187    NAVY
696.AB/BA      hP6      P6_3/mmc                 WN2_194de       #194    NAVY
697.AB/BA      tP10     P4_2/mmc                 Hagg_N4W6       #131    NAVY
698.AB/BA      oP16     Pnma                     MoO3            #62     NAVY
699.AB/BA      hP12     P6/mmm                   O3W             #191    NAVY
700.AB/BA      oP16     Pnma                     P3Tc            #62     NAVY
701.AB/BA      hP5      P-3m1                    WN4_164ad       #164    NAVY
702.AB/BA      hP20     P6_3/mmc                 B4W             #194    NAVY
703.AB/BA      oP40     Pbca                     P4Re            #61     NAVY
704.AB/BA      hP14     P6_3/mmc     B8_1v2                      #194    NAVY
705.AB/BA      oS28     Cmca         B8_1v2p                     #64     NAVY
706.AB/BA      oS14     Cmmm                     Na4Cl3v1        #65     NAVY
707.AB/BA      cP7      P-43m                    C4W3            #215    NAVY
708.AB/BA      cP7      Pm-3m                    S3U4            #221    NAVY
709.AB/BA      mS14     C2/m                     O5Nb2           #12     NAVY
710.AB/BA      hP13     P-3m1        B8_1v3                      #164    NAVY
711.AB/BA      cF52     Fm-3m                    Na7Cl6_v3       #225    NAVY
712.AB/BA      cF60     Fm-3m                    Na7Cl8v1        #225    NAVY
713.AB/BA      hP15     P-3m1        B8_1v1                      #164    NAVY
714.AB/BA      hP15     P-6m2        B8_1v1p                     #187    NAVY
715.AB/BA      tP2      P4/mmm                   Na1Cl1v2        #123    NAVY
716.AB/BA      tP6      P4/mmm                   Na3Cl3v2        #123    NAVY
717.AB/BA      hP14     P3m1         B8_1v2pp                    #156    NAVY
718.AB/BA      hR6      R3m          B3v2                        #160    NAVY
719.AB/BA      cF56     F-43m                    Na7Cl7v2        #216    NAVY
720.AB/BA      cF8      Fm-3m        B1                          #225    NAVY
721.AB/BA      hP6      P-62m                    NTa             #189    NAVY
722.AB/BA      cP6      Pm-3m                    NbO             #221    NAVY
723.AB/BA      hP16     P6_3mc                   NbS             #186    NAVY
724.AB/BA      tP3      P4/mmm                   Hagg_NW2        #123    NAVY
725.AB/BA      tI12     I4_1/amd                 Mo2N            #141    NAVY
726.AB/BA      hP9      P-31m                    NV2             #162    NAVY
727.AB/BA      hP9      P-3m1                    NW2_Pearson     #164    NAVY
728.A/B        cP8      Pm-3n        beta-W                      #223    NAVY
729.AB/BA      oP76     Pnma                     Ge8Mn11         #62     Pauling
730.AB/BA      hP48     P-3c1                    NaTe3           #165    ICSD_61355
731.AB/BA      oP48     Pbcn                     NaTe            #60     ICSD_61382
732.AB/BA      mP40     P121/c1                  Sb3Sr2          #60     ICSD_32033
733.AB/BA      mS16     C12/m1                   SrAs3           #12     ICSD_41831
734.AB/BA      tI12     I4/mmm                   La2Sb           #139    ICSD_10442
735.AB/BA      oP32     Pnma                     Ga5Tm3          #62     ICSD_104009


**** BORIDES PROJECT
1000.AB/BA     cP7      Pm-3m                    B6Si1           #221    ICSD_20240    binary_boride_ICSD_frequency=21
1001.AB/BA     hP3      P6/mmm                   B2Zr1           #191    ICSD_30327    binary_boride ICSD frequency=17
1002.AB/BA     cF52     Fm-3m                    B12Y1           #225    ICSD_23860    binary_boride_ICSD_frequency=16
1003.AB/BA     hP3      P6/mmm                   B2Mg1           #191    ICSD_26675    binary_boride_ICSD_frequency=13
1004.AB/BA     tP20     P4/mbm                   B4La1           #127    ICSD_2360     binary_boride_ICSD_frequency=12
1005.AB/BA     oS8      Cmcm                     B1Ni1           #63     ICSD_26937    binary_boride_ICSD_frequency=8
1006.AB/BA     tI12     I4/mcm                   B1Mo2           #140    ICSD_24278    binary_boride_ICSD_frequency=8
1007.AB        cF8      F-43m                    As1B1           #216    ICSD_43871    binary_boride_ICSD_frequency=5
1008.AB/BA     oI14     Immm                     B4Cr3           #71     ICSD_33704    binary_boride_ICSD_frequency=5
1009.AB        cF8      Fm-3m                    B1Pu1           #225    ICSD_43661    binary_boride_ICSD_frequency=4
1010.AB/BA     tP10     P4/mbm                   B2V3            #127    ICSD_88317    binary_boride_ICSD_frequency=4
1011.AB/BA     hP3      P6/mmm                   B2Zr1           #191    ICSD_169458   binary_boride ICSD frequency=3
1012.AB/BA     oS20     Cmcm                     B3V2            #63     ICSD_79258    binary_boride_ICSD_frequency=3
1013.AB/BA     tP20     P4/mbm                   B4Ce1           #127    ICSD_612719   binary_boride ICSD frequency=3
1014.AB/BA     cF116    Fm-3m                    B6Fe23          #225    ICSD_54786    binary_boride_ICSD_frequency=3
1015.AB        hP2      P-6m2                    B1Tc1           #187    ICSD_168895   binary_boride_ICSD_frequency=2
1016.AB/BA     hP4      P6_{3}/mmc               B1Pt1           #194    ICSD_24363    binary_boride_ICSD_frequency=2
1017.AB/BA     hP6      P6_{3}/mmc               B2Re1           #194    ICSD_23871    binary_boride_ICSD_frequency=2
1018.AB/BA     oP6      Pmmn                     B2Ru1           #59     ICSD_31871    binary_boride_ICSD_frequency=2
1019.AB/BA     mS14     C2/m                     B5W2            #12     ICSD_20326    binary_boride_ICSD_frequency=2
1020.AB/BA     oP8      Pnma                     B1Fe1           #62     ICSD_391329   binary_boride ICSD frequency=2
1021.AB/BA     tI16     I4_{1}/amd               B1Mo1           #141    ICSD_24280    binary_boride_ICSD_frequency=2
1022.AB/BA     oS16     Cmcm                     B1Re3           #63     ICSD_43662    binary_boride_ICSD_frequency=2
1023.AB/BA     hP10     P6_{3}/mmc               B3Ru2           #194    ICSD_23715    binary_boride_ICSD_frequency=2
1024.AB/BA     oP16     Pnma                     B1Ni3           #62     ICSD_24306    binary_boride_ICSD_frequency=2
1025.AB/BA     hP20     P6_{3}mc                 B3Ru7           #186    ICSD_44343    binary_boride_ICSD_frequency=2
1026.AB        cP2      Pm-3m                    B1Tc1           #221    ICSD_168900   binary_boride_ICSD_frequency=1
1027.AB/BA     cF12     Fm-3m                    B1Be2           #225    ICSD_20384    binary_boride_ICSD_frequency=1
1028.AB/BA     hP3      P6/mmm                   B2Nb1           #191    ICSD_614887   binary_boride ICSD frequency=1
1029.AB/BA     oS8      Cmcm                     B1Rh1           #63     ICSD_150732   binary_boride_ICSD_frequency=1
1030.AB/BA     hP4      P6_{3}/mmc               B1Tc1           #194    ICSD_168897   binary_boride_ICSD_frequency=1
1031.AB/BA     hP4      P6_{3}/mmc               B1Tc1           #194    ICSD_168898   binary_boride_ICSD_frequency=1
1032.AB/BA     hR4      R-3m                     B2Li2           #166    Kolmogorov-Curtarolo LiB    binary_boride_ICSD_frequency=1
1033.AB/BA     oI10     Immm                     B4Cr1           #71     ICSD_24353    binary_boride ICSD frequency=1
1034.AB/BA     mS10     C2/m                     B4Mn1           #12     ICSD_15079    binary_boride_ICSD_frequency=1
1035.AB/BA     tI12     I4/mcm                   B1Fe2           #140    ICSD_160788   binary_boride ICSD frequency=1
1036.AB/BA     tI12     I4/mcm                   B1Fe2           #140    ICSD_160789   binary_boride_ICSD_frequency=1
1037.AB/BA     tI12     I-42m                    B1Fe2           #121    ICSD_16809    binary_boride_ICSD_frequency=1
1038.AB/BA     oP6      Pnnm                     B1Pd2           #58     ICSD_10487    binary_boride_ICSD_frequency=1
1039.AB/BA     hP6      P6_{3}/mmc               B1Pt2           #194    ICSD_615207   binary_boride_ICSD_frequency=1
1040.AB/BA     tI12     I4/mcm                   B1Ti2           #140    ICSD_189385   binary_boride_ICSD_frequency=1
1041.AB/BA     mS12     C2/m                     B2Mo1           #12     ICSD_418398   binary_boride ICSD frequency=1
1042.AB/BA     cF24     Fd-3m                    B2Sc1           #227    ICSD_106903   binary_boride_ICSD_frequency=1
1043.AB/BA     hR6      R-3m                     B4Mo2           #166    ICSD_39554    binary_boride_ICSD_frequency=1
1044.AB/BA     oI14     Immm                     B4Mn3           #71     ICSD_614732   binary_boride ICSD frequency=1
1045.AB/BA     oI14     Immm                     B4Nb3           #71     ICSD_76631    binary_boride ICSD frequency=1
1046.AB/BA     mS14     C2/m                     B5W2            #12     ICSD_20326    binary_boride_ICSD_frequency=1
1047.AB/BA     oP8      Pnma                     B1Co1           #62     ICSD_612863   binary_boride_ICSD_frequency=1
1048.AB/BA     hP8      P6_{3}22                 B1Fe3           #182    ICSD_184958   binary_boride_ICSD_frequency=1
1049.AB/BA     oP8      Pnma                     B1Li1           #62     ICSD_153291   binary_boride_ICSD_frequency=1
1050.AB/BA     oP8      Pnma                     B1Mn1           #62     ICSD_76630    binary_boride ICSD frequency=1
1051.AB/BA     hP8      P6_{3}/mmc               B2Li2           #194    Kolmogorov-Curtarolo LiB    binary_boride_ICSD_frequency=1
1052.AB/BA     hR8      R-3m                     B3Mo1           #166    ICSD_167734   binary_boride_ICSD_frequency=1
1053.AB/BA     hP8      P6_{3}/mmc               B3Re1           #194    ICSD_24361    binary_boride_ICSD_frequency=1
1054.AB/BA     hP10     P6_{3}/mmc               B3Ru2           #194    ICSD_108082   binary_boride_ICSD_frequency=1
1055.AB/BA     oP10     Pnnm                     B4Cr1           #58     ICSD_186851   binary_boride_ICSD_frequency=1
1056.AB/BA     oP10     Pnnm                     B4Fe1           #58     ICSD_425311   binary_boride ICSD frequency=1
1057.AB/BA     hP10     P6_{3}/mmc               B4Mo1           #194    ICSD_167735   binary_boride_ICSD_frequency=1
1058.AB/BA     oS22     Cmmm                     B6Ta5           #65     ICSD_68538    binary_boride_ICSD_frequency=1
1059.AB/BA     cF48     F-43m                    B11Li1          #216    ICSD_164841   binary_boride_ICSD_frequency=1
1060.AB/BA     hR12     R3m                      B11Li1          #160    ICSD_164842   binary_boride_ICSD_frequency=1
1061.AB/BA     cF48     F-43m                    B11Li1          #216    ICSD_164843   binary_boride_ICSD_frequency=1
1062.AB/BA     tI24     I-4m2                    B11Li1          #119    ICSD_164844   binary_boride_ICSD_frequency=1
1063.AB/BA     tI24     I-4m2                    B11Li1          #119    ICSD_164845   binary_boride_ICSD_frequency=1
1064.AB/BA     oF48     Fddd                     B1Cr2           #70     ICSD_188474   binary_boride_ICSD_frequency=1
1065.AB/BA     oP12     Pnma                     B1Rh2           #62     ICSD_24699    binary_boride_ICSD_frequency=1
1066.AB/BA     oP12     Pnma                     B2Fe1           #62     ICSD_425310   binary_boride_ICSD_frequency=1
1067.AB/BA     hP12     P6_{3}/mmc               B2W1            #194    ICSD_23716    binary_boride_ICSD_frequency=1
1068.AB/BA     tI26     I4/mmm                   B12Sc1          #139    ICSD_615424   binary_boride_ICSD_frequency=1
1069.AB/BA     mS28     C2/c                     B2Pd5           #15     ICSD_43513    binary_boride_ICSD_frequency=1
1070.AB/BA     mS28     C2/c                     B3Ni4           #15     ICSD_150561   binary_boride ICSD frequency=1
1071.AB/BA     mS28     C2/c                     B3Ni4           #15     ICSD_24308    binary_boride_ICSD_frequency=1
1072.AB/BA     hP14     P6_{3}/mmc               B5Mo2           #194    ICSD_167733   binary_boride_ICSD_frequency=1
1073.AB/BA     hP14     P6_{3}/mmc               B5W2            #194    ICSD_615696   binary_boride_ICSD_frequency=1
1074.AB/BA     hR15     R-3m                     B12Si3          #166    ICSD_615435   binary_boride ICSD frequency=1
1075.AB/BA     oP16     Pnma                     B1Co3           #62     ICSD_44339    binary_boride_ICSD_frequency=2
1076.AB/BA     oP16     Pnma                     B1Pd3           #62     ICSD_43514    binary_boride_ICSD_frequency=1
1077.AB/BA     tP16     P4_{2}/mmc               B4Mn4           #131    ICSD_74985    binary_boride_ICSD_frequency=1
1078.AB/BA     hP18     R-3m                     B2Mo1           #166    ICSD_418397   binary_boride ICSD frequency=1
1079.AB/BA     hP18     P6_{3}/mmc               B4Rh5           #194    ICSD_86395    binary_boride_ICSD_frequency=1
1080.AB/BA     hP20     P6_{3}/mmc               B4W1            #194    ICSD_43193    binary_boride_ICSD_frequency=1
1081.AB/BA     tP20     P4/mbm                   B4W1            #127    ICSD_615683   binary_boride_ICSD_frequency=1
1082.AB/BA     hP20     P6_{3}cm                 B9Li1           #185    ICSD_164625   binary_boride_ICSD_frequency=1
1083.AB/BA     oS46     Cmmm                     B20Na3          #65     ICSD_407159   binary_boride_ICSD_frequency=1
1084.AB/BA     oP28     Pnma                     B3Ni4           #62     ICSD_24307    binary_boride_ICSD_frequency=1
1085.AB/BA     oI64     Imma                     B3Si1           #74     ICSD_412621   binary_boride_ICSD_frequency=1
1086.AB/BA     oP38     Pbam                     B8Ru11          #55     ICSD_43663    binary_boride_ICSD_frequency=1


**** CARBIDE PROJECT
001.AA         fcc/cF4  Fm-3m         A1         Cu              #225    from gndstate
001.BB         fcc/cF4  Fm-3m         A1         Cu              #225    from gndstate
058.AA         bcc/cI2  Im-3m         A2         W               #229    from gndstate
058.BB         bcc/cI2  Im-3m         A2         W               #229    from gndstate
115.AA         hcp/hcp  P6_3/mmc      A3         Mg              #194    from gndstate
115.BB         hcp/hcp  P6_3/mmc      A3         Mg              #194    from gndstate
201.AB         cF8      Fm-3m         B1         NaCl            #225    from gndstate
061.AB         cP2      Pm-3m         B2         CsCl            #221    from gndstate
218.AB         cF8      F-43m         B3         ZnS             #216    from gndstate
219.AB         hP4      P63/mc        B4         ZnS             #186    from gndstate
116.AB         hP2      P_6m2         B_h        WC              #187    from gndstate
301.AA         cF8      Fd-3m         A4         C (DIA)         #227    BOOK
301.BB         cF8      Fd-3m         A4         C (DIA)         #227    BOOK
325.AA         hP4      P63/mmc       A9         C (GRA)         #194    BOOK
325.BB         hP4      P63/mmc       A9         C (GRA)         #194    BOOK
1200.AB        cF8      Fm-3m                    C1Ti1           #225    ICSD_1546     binary_carbide_ICSD_frequency=30
1201.AB        cF8      F-43m                    C1Si1           #216    ICSD_24217    binary_carbide_ICSD_frequency=20
1202.AB/BA     tI6      I4/mmm                   C2U1            #139    ICSD_14330    binary_carbide_ICSD_frequency=18
1203.AB/BA     cI40     I-43d                    C3Pr2           #220    ICSD_2444     binary_carbide_ICSD_frequency=17
1204.AB        hP2      P-6m2                    C1W1            #187    ICSD_15406    binary_carbide_ICSD_frequency=12
1205.AB        cP2      Pm-3m                    C1Sc1           #221    ICSD_181048   binary_carbide_ICSD_frequency=11
1206.AB/BA     cI132    Im-3                     C60K6           #204    ICSD_66879    binary_carbide_ICSD_frequency=5
1207.AB/BA     hR3      R-3m                     C1Y2            #166    ICSD_22283    binary_carbide_ICSD_frequency=4
1208.AB/BA     mS12     C2/c                     C2Ca1           #15     ICSD_54184    binary_carbide_ICSD_frequency=4
1209.AB/BA     oP16     Pnma                     C1Fe3           #62     ICSD_16593    binary_carbide_ICSD_frequency=4
1210.AB/BA     tP70     P4/mnc                   C4Sc3           #128    ICSD_71145    binary_carbide_ICSD_frequency=4
1211.AB/BA     tI6      I4/mmm                   Ba1C2           #139    ICSD_56160    binary_carbide_ICSD_frequency=3
1212.AB/BA     hP3      P-3m1                    C1Nb2           #164    ICSD_43670    binary_carbide_ICSD_frequency=3
1213.AB/BA     oP12     Pbcn                     C1Mo2           #60     ICSD_65701    binary_carbide_ICSD_frequency=3
1214.AB/BA     cP12     Pa-3                     C2Th1           #205    ICSD_77470    binary_carbide_ICSD_frequency=3
1215.AB/BA     oP16     Pnma                     C2Rb2           #62     ICSD_51529    binary_carbide_ICSD_frequency=3 
1216.AB/BA     oP18     Pbam                     C5Y4            #55     ICSD_84869    binary_carbide_ICSD_frequency=3
1217.AB/BA     tP68     P-42_{1}c                C19Sc15         #114    ICSD_42631    binary_carbide_ICSD_frequency=1
1218.AB/BA     cF12     Fm-3m                    Be2C1           #225    ICSD_41567    binary_carbide_ICSD_frequency=2
1219.AB/BA     hP3      P-6m2                    C2Os1           #187    ICSD_168279   binary_carbide_ICSD_frequency=2
1220.AB/BA     oI8      Immm                     C2Li2           #71     ICSD_25705    binary_carbide_ICSD_frequency=2
1221.AB/BA     hP6      P6_{3}/mmc               C2Os1           #194    ICSD_168280   binary_carbide_ICSD_frequency=2
1222.AB/BA     oP12     Pbcn                     C1V2            #60     ICSD_9965     binary_carbide_ICSD_frequency=2
1223.AB/BA     mP22     P2_{1}/c                 C7Ho4           #14     ICSD_83381    binary_carbide_ICSD_frequency=2
1224.AB/BA     cF116    Fm-3m                    C6Cr23          #225    ICSD_69536    binary_carbide_ICSD_frequency=2
1225.AB/BA     hR3      R-3m                     Ba1C2           #166    ICSD_186576   binary_carbide_ICSD_frequency=1
1226.AB/BA     hP3      P-6m2                    C1Mo2           #187    ICSD_43669    binary_carbide_ICSD_frequency=1
1227.AB/BA     hP3      P-3m1                    C1Sc2           #164    ICSD_280743   binary_carbide_ICSD_frequency=1
1228.AB/BA     tI6      I4/mmm                   C2Re1           #139    ICSD_184661   binary_carbide_ICSD_frequency=1
1229.AB/BA     tP3      P4/mmm                   C2Re1           #123    ICSD_184664   binary_carbide_ICSD_frequency=1
1230.AB        hP4      P6_{3}mc                 C1Nb1           #186    ICSD_189092   binary_carbide_ICSD_frequency=1
1231.AB/BA     hP4      P6_{3}/mmc               C1Nb1           #194    ICSD_189093   binary_carbide_ICSD_frequency=1
1232.AB/BA     hP4      P6_{3}/mmc               C1Os1           #194    ICSD_168278   binary_carbide_ICSD_frequency=1
1233.AB/BA     hR4      R-3m                     C1Ru1           #166    ICSD_188285   binary_carbide_ICSD_frequency=1
1234.AB/BA     hP4      P6_{3}/mmc               C1Sc1           #194    ICSD_189088   binary_carbide_ICSD_frequency=1
1235.AB        hP4      P6_{3}mc                 C1Si1           #186    ICSD_24261    binary_carbide_ICSD_frequency=1
1236.AB/BA     tI10     I4/m                     C1Fe4           #87     ICSD_187144   binary_carbide_ICSD_frequency=1
1237.AB/BA     cP5      P-43m                    C1Fe4           #215    ICSD_44729    binary_carbide_ICSD_frequency=1
1238.AB/BA     oI12     Imm2                     B1C5            #44     ICSD_166554   binary_carbide_ICSD_frequency=1
1239.AB/BA     hP6      P3m1                     B1C5            #156    ICSD_166555   binary_carbide_ICSD_frequency=1
1240.AB/BA     oP6      Pnnm                     C1Co2           #58     ICSD_617391   binary_carbide_ICSD_frequency=1
1241.AB/BA     hP3      P6/mmm                   C1Fe2           #191    ICSD_162103   binary_carbide_ICSD_frequency=1
1242.AB/BA     oP6      Pnnm                     C1Fe2           #58     ICSD_76826    binary_carbide_ICSD_frequency=2
1243.AB/BA     hP6      P6_{3}/mmc               C1Re2           #194    ICSD_180376   binary_carbide_ICSD_frequency=1
1244.AB/BA     mS12     C2/m                     C2Ca1           #12     ICSD_54185    binary_carbide_ICSD_frequency=1
1245.AB/BA     mS12     C2/m                     C2Ca1           #12     ICSD_54188    binary_carbide_ICSD_frequency=1
1246.AB/BA     mS12     C2/m                     C2Ir1           #12     ICSD_181488   binary_carbide_ICSD_frequency=1
1247.AB/BA     tI12     I4/mcm                   C2Ir1           #140    ICSD_181489   binary_carbide_ICSD_frequency=1
1248.AB/BA     tP6      P4_{2}/mnm               C2Mg1           #136    ICSD_88057    binary_carbide_ICSD_frequency=1
1249.AB/BA     hP6      P6_{3}/mmc               C2Re1           #194    ICSD_184660   binary_carbide_ICSD_frequency=1
1250.AB/BA     tP6      P4_{2}/mmc               C2Si1           #131    ICSD_187721   binary_carbide_ICSD_frequency=1
1251.AB/BA     mS12     C2/c                     C2Th1           #15     ICSD_202306   binary_carbide_ICSD_frequency=1
1252.AB/BA     mS12     C2/c                     C2Th1           #15     ICSD_26568    binary_carbide_ICSD_frequency=1
1253.AB/BA     tP6      P4_{2}/mmc               C2Th1           #131    ICSD_77469    binary_carbide_ICSD_frequency=1
1254.AB/BA     hR7      R3m                      Al4C3           #160    ICSD_14397    binary_carbide_ICSD_frequency=1
1255.AB/BA     cP7      P-43m                    As4C3           #215    ICSD_163831   binary_carbide_ICSD_frequency=1
1256.AB/BA     cP7      Pm-3m                    C3Nb4           #221    ICSD_42758    binary_carbide_ICSD_frequency=1
1257.AB/BA     tP8      P-4m2                    B1C7            #115    ICSD_181952   binary_carbide_ICSD_frequency=1
1258.AB/BA     hP8      P3m1                     B1C7            #156    ICSD_181953   binary_carbide_ICSD_frequency=1
1259.AB/BA     cP8      P-43m                    B1C7            #215    ICSD_181954   binary_carbide_ICSD_frequency=1
1260.AB/BA     oP8      Pmm2                     B1C7            #25     ICSD_181955   binary_carbide_ICSD_frequency=1
1261.AB/BA     hR8      R3m                      B1C7            #160    ICSD_181956   binary_carbide_ICSD_frequency=1
1262.AB/BA     hP8      P6_{3}22                 C1Fe3           #182    ICSD_163151   binary_carbide_ICSD_frequency=1
1263.AB/BA     hP8      P6_{3}22                 C1Fe3           #182    ICSD_42542    binary_carbide_ICSD_frequency=1
1264.AB/BA     hP8      P6_{3}/mmc               C1Mo1           #194    ICSD_618298   binary_carbide_ICSD_frequency=1
1265.AB/BA     hR8      R-3c                     C1Ni3           #167    ICSD_17005    binary_carbide_ICSD_frequency=1
1266.AB/BA     cP8      P2_{1}3                  C1Os1           #198    ICSD_168277   binary_carbide_ICSD_frequency=1
1267.AB/BA     hP8      P6_{3}/mmc               C1Re1           #194    ICSD_618702   binary_carbide_ICSD_frequency=1
1268.AB        hP8      P6_{3}mc                 C1Si1           #186    ICSD_24170    binary_carbide_ICSD_frequency=1
1269.AB/BA     hP9      P-31m                    C1W2            #162    ICSD_162419   binary_carbide_ICSD_frequency=1
1270.AB/BA     hP9      P6/mmm                   C8Cs1           #191    ICSD_74641    binary_carbide_ICSD_frequency=1
1271.AB/BA     mP10     P2_{1}/m                 C1Fe4           #11     ICSD_187143   binary_carbide_ICSD_frequency=1
1272.AB        hP10	    P3m1                     C1Si1           #156    ICSD_107204   binary_carbide_ICSD_frequency=1
1273.AB        hR10    	R3m                      C1Si1           #160    ICSD_24168    binary_carbide_ICSD_frequency=1
1274.AB/BA     oS20     Cmcm                     C2Cr3           #63     ICSD_57009    binary_carbide_ICSD_frequency=1
1275.AB/BA     oP10     Pnnm                     C3Mg2           #58     ICSD_71941    binary_carbide_ICSD_frequency=1
1276.AB/BA     mP10     P2_{1}/c                 C4Ir1           #14     ICSD_181493   binary_carbide_ICSD_frequency=1
1277.AB/BA     oS20     Ccce                     C4Ir1           #68     ICSD_181496   binary_carbide_ICSD_frequency=1
1278.AB/BA     tP10     P4/nbm                   C4Ir1           #125    ICSD_181497   binary_carbide_ICSD_frequency=1
1279.AB/BA     mS22     C2/m                     C5Nb6           #12     ICSD_20695    binary_carbide_ICSD_frequency=1
1280.AB/BA     tI12     I-4m2                    B1C5            #119    ICSD_166553   binary_carbide_ICSD_frequency=1
1281.AB/BA     aP12     P-1                      B1C5            #2      ICSD_166556   binary_carbide_ICSD_frequency=1
1282.AB/BA     oP12     Pnma                     C1Fe2           #62     ICSD_187138   binary_carbide_ICSD_frequency=1
1283.AB/BA     hP12     P6_{3}/mmc               C1Mo1           #194    ICSD_44987    binary_carbide_ICSD_frequency=1
1284.AB/BA     oP12     Pbcn                     C1Mo2           #60     ICSD_43322    binary_carbide_ICSD_frequency=1
1285.AB/BA     oP12     Pnma                     C1Nb2           #62     ICSD_31973    binary_carbide_ICSD_frequency=1
1286.AB/BA     oP12     Pmma                     C1Nb2           #51     ICSD_653912   binary_carbide_ICSD_frequency=1
1287.AB        hP12 	  P6_{3}mc                 C1Si1           #186    ICSD_15325    binary_carbide_ICSD_frequency=1
1288.AB/BA     cF48     Fd-3m                    C1Ti2           #227    ICSD_77473    binary_carbide_ICSD_frequency=1
1289.AB/BA     mP12     P2_{1}/c                 C2Ca1           #14     ICSD_24074    binary_carbide_ICSD_frequency=1
1290.AB/BA     cP12     Pn-3m                    C2Ca1           #224    ICSD_31092    binary_carbide_ICSD_frequency=1
1291.AB/BA     aP12     P-1                      C2Ca1           #2      ICSD_66663    binary_carbide_ICSD_frequency=1
1292.AB/BA     hP12     P-62m                    C2Cs2           #189    ICSD_51533    binary_carbide_ICSD_frequency=1
1293.AB/BA     cP12     P4_{3}32                 C2Na1           #212    ICSD_186177   binary_carbide_ICSD_frequency=1
1294.AB/BA     hR13     R-3m                     C5Ti8           #166    ICSD_20822    binary_carbide_ICSD_frequency=1
1295.AB        hR14     R3m                      C1Si1           #160    ICSD_42513    binary_carbide_ICSD_frequency=1
1296.AB/BA     mS28     C2/c                     C2Fe5           #15     ICSD_423885   binary_carbide_ICSD_frequency=1
1297.AB/BA     mS28     C2/c                     C2Fe5           #15     ICSD_43194    binary_carbide_ICSD_frequency=1
1298.AB/BA     mS28     C2/c                     C2Fe5           #15     ICSD_89328    binary_carbide_ICSD_frequency=1
1299.AB/BA     mS28     C2/c                     C2Mn5           #15     ICSD_618246   binary_carbide_ICSD_frequency=1
1300.AB/BA     mS28     C2/c                     C2Mn5           #15     ICSD_69535    binary_carbide_ICSD_frequency=1
1301.AB/BA     cI28     I-43d                    C3Sc4           #220    ICSD_42760    binary_carbide_ICSD_frequency=1
1302.AB/BA     hP14     P6_{3}/mmc               C6Eu1           #194    ICSD_169041   binary_carbide_ICSD_frequency=1
1303.AB/BA     hP14     P6_{3}/mmc               C6Yb1           #194    ICSD_601565   binary_carbide_ICSD_frequency=1
1304.AB/BA     hR15     R-3m                     B13C2           #166    ICSD_446      binary_carbide_ICSD_frequency=1
1305.AB/BA     hR15     R-3m                     B4C1            #166    ICSD_29093    binary_carbide_ICSD_frequency=1
1306.AB/BA     oP16     Pnma                     C1Cr3           #62     ICSD_617486   binary_carbide_ICSD_frequency=1
1307.AB/BA     oP16     Pnma                     C1Fe3           #62     ICSD_167667   binary_carbide_ICSD_frequency=1
1308.AB/BA     tI32     I-4                      C1Fe3           #82     ICSD_187142   binary_carbide_ICSD_frequency=1
1309.AB/BA     oP16     Pnma                     C1Fe3           #62     ICSD_99010    binary_carbide_ICSD_frequency=1
1310.AB        hR16    	R3m                      C1Si1           #160    ICSD_18136    binary_carbide_ICSD_frequency=1
1311.AB        hP16     P6_{3}mc                 C1Si1           #186    ICSD_24630    binary_carbide_ICSD_frequency=1
1312.AB/BA     cF64     Fd-3m                    C1Ti1           #227    ICSD_618950   binary_carbide_ICSD_frequency=1
1313.AB/BA     oP16     Pnma                     C2Cs2           #62     ICSD_51531    binary_carbide_ICSD_frequency=1
1314.AB/BA     tI32     I4_{1}/acd               C2K2            #142    ICSD_36142    binary_carbide_ICSD_frequency=1
1315.AB/BA     tI32     I4_{1}/acd               C2Na2           #142    ICSD_28066    binary_carbide_ICSD_frequency=1
1316.AB/BA     oP16     Pnma                     C3Ir1           #62     ICSD_181492   binary_carbide_ICSD_frequency=1
1317.AB/BA     oP16     Pbam                     C3Ir5           #55     ICSD_181485   binary_carbide_ICSD_frequency=1
1318.AB        hR18     R3m                      C1Si1           #160    ICSD_24631    binary_carbide_ICSD_frequency=1
1319.AB/BA     oF72     Fddd                     C8K1            #70     ICSD_70020    binary_carbide_ICSD_frequency=1
1320.AB/BA     oF72     Fddd                     C8Rb1           #70     ICSD_200563   binary_carbide_ICSD_frequency=1
1321.AB        hP20     P3m1                     C1Si1           #156    ICSD_38150    binary_carbide_ICSD_frequency=1
1322.AB/BA     oP20     Pnma                     C2Cr3           #62     ICSD_15086    binary_carbide_ICSD_frequency=1
1323.AB/BA     oP20     Pnma                     C2Cr3           #62     ICSD_653958   binary_carbide_ICSD_frequency=1
1324.AB/BA     hP20     P6_{3}mc                 C3Fe7           #186    ICSD_76830    binary_carbide_ICSD_frequency=1
1325.AB/BA     tP20     P-4c2                    C3Os2           #116    ICSD_168275   binary_carbide_ICSD_frequency=1
1326.AB/BA     oS40     P6_{3}/mmc               C4Ir1           #194    ICSD_181498   binary_carbide_ICSD_frequency=1
1327.AB/BA     cI40     Im-3m                    C7Ir3           #229    ICSD_181491   binary_carbide_ICSD_frequency=1
1328.AB/BA     mP22     P2_{1}/c                 C7Lu4           #14     ICSD_83382    binary_carbide_ICSD_frequency=1
1329.AB/BA     oP24     Pnma                     Ba1C2           #62     ICSD_168408   binary_carbide_ICSD_frequency=1
1330.AB/BA     cF112    Fd-3m                    C1Mn6           #227    ICSD_187037   binary_carbide_ICSD_frequency=1
1331.AB        hP28     P3m1                     C1Si1           #156    ICSD_28303    binary_carbide_ICSD_frequency=1
1332.AB/BA     cF116    Fm-3m                    C6Mn23          #225    ICSD_154720   binary_carbide_ICSD_frequency=1
1333.AB        hP32     P6_{3}mc                 C1Si1           #186    ICSD_42858    binary_carbide_ICSD_frequency=1
1334.AB/BA     tI64     I4_{1}/acd               C2K2            #142    ICSD_89528    binary_carbide_ICSD_frequency=1
1335.AB/BA     hP33     P3_{1}12                 C5V6            #151    ICSD_71098    binary_carbide_ICSD_frequency=1
1336.AB        hP36     P3m1                     C1Si1           #156    ICSD_42859    binary_carbide_ICSD_frequency=1
1337.AB/BA     tP36     P-4c2                    C5Ir4           #116    ICSD_181486   binary_carbide_ICSD_frequency=1
1338.AB        hP38     P3m1                     C1Si1           #156    ICSD_27635    binary_carbide_ICSD_frequency=1
1339.AB/BA     oP40     Pbcn                     C3Os2           #60     ICSD_168276   binary_carbide_ICSD_frequency=1
1340.AB        hP42     P3m1                     C1Si1           #156    ICSD_43827    binary_carbide_ICSD_frequency=1
1341.AB/BA     mP44     P2_{1}/c                 C7Ho4           #14     ICSD_74817    binary_carbide_ICSD_frequency=1
1342.AB/BA     cP60     P4_{3}32                 C0.88V1         #212    ICSD_60191    binary_carbide_ICSD_frequency=1
1343.AB/BA     cP60     P4_{3}32                 C7V8            #212    ICSD_22177    binary_carbide_ICSD_frequency=1
1344.AB/BA     cP60     P4_{3}32                 C7V8            #212    ICSD_43259    binary_carbide_ICSD_frequency=1
1345.AB/BA     cP60     P4_{3}32                 C7V8            #212    ICSD_85953    binary_carbide_ICSD_frequency=1
1346.AB/BA     tP68     P-42_{1}c                C19Er15         #114    ICSD_76817    binary_carbide_ICSD_frequency=3
1347.AB/BA     cP126    Pm-3n                    Ba3C60          #223    ICSD_70063    binary_carbide_ICSD_frequency=1
1348.AB/BA     cP126    Pm-3                     C60Sr3          #200    ICSD_75354    binary_carbide_ICSD_frequency=1


*** TERNARIES
T0001.A2BC/AB2C/ABC2              cF16  Fm(-3)m  AlCu2Mn-Heusler              #225    Navy [Wyckoff A~C, 3 protos] [all --aflow_proto T0001 A B C]
T0001[-tet[+,-]*%].A2BC/AB2C/ABC2 cF16  Fm(-3)m  AlCu2Mn-Heusler-tet*%        #225    Navy [Wyckoff A~C, 3 protos] [all --aflow_proto T0001 A B C] TET distorted along z +-*%
T0002.A2BC/AB2C/ABC2              cF16  F-43m    Cu1Li2Sn1-antiHeusler        #216    ICSD_15130 [Wyckoff A~C, 3 protos]
T0002[-tet[+,-]*%].A2BC/AB2C/ABC2 cF16  F-43m    Cu1Li2Sn1-antiHeusler-tet*%  #216    ICSD_15130 [Wyckoff A~C, 3 protos]
T0001.{ABC}                       cF16  Fm(-3)m  AlCu2Mn-Heusler              #225    Navy [Wyckoff A~C, 3 protos] [all --aflow_proto T0001 A B C ]
T0001[-tet[+,-]*%].{ABC}          cF16  Fm(-3)m  AlCu2Mn-Heusler-tet*%        #225    Navy [Wyckoff A~C, 3 protos] [all --aflow_proto T0001 A B C ] TET distorted along Z
T0002.{ABC}                       cF16  F-43m    Cu1Li2Sn1-antiHeusler        #216    ICSD_15130 [Wyckoff A~C, 3 protos]
T0002[-tet[+,-]*%].{ABC}          cF16  F-43m    Cu1Li2Sn1-antiHeusler-tet*%  #216    ICSD_15130 [Wyckoff A~C, 3 protos]  TET distorted along Z
T0003.{ABC}                       cF12  F-43m    AgAsMg-halfHeusler       C1b #216    ICSD_43819
T0003[-tet[+,-]*%].{ABC}          cF12  F-43m    AgAsMg-halfHeusler-tet*% C1b #216    ICSD_43819  TET distorted along Z
T0004.{ABC}    oP12     Pnma            Co1Ge1Mn1_ICSD_52968     #62     ICSD_52968
T0005.{ABC}    hP6      P63/mmc         Co1Ge1Mn1_ICSD_623495    #194    ICSD_623495
T0006.{ABC}    hP12     P63/mmc         Co3Ge1Mn2_ICSD_52972     #194    ICSD_52972
T0007.{ABC}    oP12     Pnma            Cl1La1Se1_ICSD_425686    #62     ICSD_425686
T0008.{ABC}    oP12     Pnma          Cd1Na1P1 Co2Si1_ICSD_44858 #62     ICSD_44858
T0009.{ABC}    cP5      Pm-3m                    perovskite      #221 
T0010.{ABC}    cF24     F-43m           Cu4MgSn       C15b       #216    ICSD_108141
T0011.{ABC}    tI10    I4/mmm           As2Ce1Pd2                #139    ICSD_604354

*** TERNARY SQS
TS001.ABC/BCA/CBA  oS64    Cmma                 SQS_L12          #67     SQS_L12

*** TERNARY FCC
TFCC001.ABC           hP3       P3m1         ternary_fcc_1        #156    3cell
TFCC002.ABC           oI6       Imm2         ternary_fcc_2        #44     3cell
TFCC003.ABC           tI6       I4mm         ternary_fcc_3        #107    3cell
TFCC004.ABC/BCA/CAB   hR4       R3m          ternary_fcc_4        #160    4cell
TFCC005.ABC/BCA/CAB   oP4       Pmm2         ternary_fcc_5        #25     4cell
TFCC006.ABC/BCA/CAB   mS8       Cm           ternary_fcc_6        #8      4cell
TFCC007.ABC/BCA/CAB   tP4       P4mm         ternary_fcc_7        #99     4cell
TFCC008.ABC/BCA/CAB   tI8       I-4m2        ternary_fcc_8        #119    4cell
TFCC009.ABC/BCA/CAB   hR4       R-3m         ternary_fcc_9        #166    4cell
TFCC010.ABC/BCA/CAB   oP4       Pmmm         ternary_fcc_10       #47     4cell
TFCC011.ABC/BCA/CAB   mS8       C2/m         ternary_fcc_11       #12     4cell
TFCC012.ABC/BCA/CAB   tP4       P4/mmm       ternary_fcc_12       #123    4cell
TFCC013.ABC/BCA/CAB   tI8       I4/mmm       ternary_fcc_13       #139    4cell
TFCC014.ABC/BCA/CAB   mP4       P2/m         ternary_fcc_14       #10     4cell
TFCC015.ABC/BCA/CAB   oS8       Cmmm         ternary_fcc_15       #65     4cell
TFCC016.ABC/BCA/CAB   tP4       P4/mmm       ternary_fcc_16       #123    4cell

*** TERNARY BCC
TBCC001.ABC           oF12      Fmm2         ternary_bcc_1        #42     3cell
TBCC002.ABC           hP3       P3m1         ternary_bcc_2        #156    3cell
TBCC003.ABC           tI6       I4mm         ternary_bcc_3        #107    3cell
TBCC004.ABC/BCA/CAB   oS8       Cmm2         ternary_bcc_4        #35     4cell
TBCC005.ABC/BCA/CAB   mP4       Pm           ternary_bcc_5        #6      4cell
TBCC006.ABC/BCA/CAB   tP4       P4mm         ternary_bcc_6        #99     4cell
TBCC007.ABC/BCA/CAB   oI8       Imm2         ternary_bcc_7        #44     4cell
TBCC008.ABC/BCA/CAB   cF16      F-43m        ternary_bcc_8        #216    4cell same as T0002 (anti-Heusler)
TBCC009.ABC/BCA/CAB   oS8       Cmmm         ternary_bcc_9        #65     4cell
TBCC010.ABC/BCA/CAB   mP4       P2/m         ternary_bcc_10       #10     4cell
TBCC011.ABC/BCA/CAB   tP4       P4/mmm       ternary_bcc_11       #123    4cell
TBCC012.ABC/BCA/CAB   oI8       Cmmm         ternary_bcc_12       #65     4cell
TBCC013.ABC/BCA/CAB   cF16      Fm-3m        ternary_bcc_13       #225    4cell same as T0001 (Heusler)
TBCC014.ABC/BCA/CAB   mS8       C2/m         ternary_bcc_14       #12     4cell
TBCC015.ABC/BCA/CAB   oP4       Pmmm         ternary_bcc_15       #47     4cell
TBCC016.ABC/BCA/CAB   tP4       P4/mmm       ternary_bcc_16       #123    4cell

*** QUATERNARY
Q0001.{ABCD}   cF40     Fm(-3)m         Al1F6K2Na1_ICSD_6027     #225    ICSD_6027 Elpasolite [Wyckoff all distinct, 24 protos]

*****************************************************************************************************************************************************
*****************************************************************************************************************************************************
ANRL PROTOTYPES
*****************************************************************************************************************************************************
ANRL - UNARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes (part)
----------------------------------------------------------------------------------------------------------------------------------------------------
A_aP4_2_aci-001                                  aP4       P-1                 -         Cf                  #2        Cf (part 1)         
A_mP4_4_2a-001                                   mP4       P2_{1}              -         Te                  #4        High-Pressure Te (part 1)
A_mC12_5_3c-001                                  mC12      C2                  A19       Po                  #5        Po (part 1)         
A_mP8_10_2m2n-001                                mP8       P2/m                -         C                   #10       S-carbon (part 2)   
A_mP16_11_8e-001                                 mP16      P2_{1}/m            -         alpha-Pu            #11       alpha-Pu (part 1)    
A_mC16_12_4i-001                                 mC16      C2/m                -         C                   #12       M-carbon (part 2)   
A_mC34_12_ah3i2j-001                             mC34      C2/m                -         beta-Pu             #12       beta-Pu (part 1)     
A_mC4_12_i-001                                   mC4       C2/m                -         O2                  #12       alpha-O2 (part 1)    
A_mP84_13_21g-001                                mP84      P2/c                -         P                   #13       Monoclinic Phosphorus (part 1)
A_mP32_14_8e-001                                 mP32      P2_{1}/c            A_l       Se                  #14       betaSe (part 1)     
A_mP64_14_16e-001                                mP64      P2_{1}/c            A_k       Se                  #14       Se (part 1)         
A_mC24_15_2e2f-001                               mC24      C2/c                -         H                   #15       H-III (part 2)      
A_oP16_55_2g2h-001                               oP16      Pbam                -         C                   #55       R-carbon (part 2)   
A_oP8_62_2c-001                                  oP8       Pnma                A_c       alpha-Np            #62       alpha-Np (part 1)    
A_oC4_63_c-001                                   oC4       Cmcm                A20       alpha-U             #63       alpha-Uranium (part 1)
A_oC8_64_f-001                                   oC8       Cmce                A11       alpha-Ga            #64       alpha-Gallium (part 1)
A_oC8_64_f-002                                   oC8       Cmce                A17       P                   #64       Black Phosphorus (part 1)
A_oC8_64_f-003                                   oC8       Cmce                A14       I2                  #64       Molecular Iodine (part 1)
A_oF128_70_4h-001                                oF128     Fddd                A16       S                   #70       alpha-Sulfur (part 1)
A_oF8_70_a-001                                   oF8       Fddd                -         gamma-Pu            #70       gamma-Pu (part 1)   
A_tP12_96_ab-001                                 tP12      P4_{3}2_{1}2        -         Si                  #96       "ST12" of Si (part 1)
A_tP4_129_ac-001                                 tP4       P4/nmm              A_d       beta-Np             #129      beta-Np (part 1)     
A_tP50_134_b2m2n-001                             tP50      P4_{2}/nnm          A_g       B                   #134      T-50 Boron (part 1) 
A_tP30_136_bf2ij-001                             tP30      P4_{2}/mnm          A_b       beta-U              #136      beta-U (part 1)      
A_tP4_136_f-001                                  tP4       P4_{2}/mnm          -         N2                  #136      gamma-N2 (part 1)    
sigma_tP30_136_bf2ij-001                         tP30      P4_{2}/mnm          D8_b      CrFe                #136      sigma-CrFe (part 1)  
A_tP12_138_bi-001                                tP12      P4_{2}/ncm          -         C                   #138      C (part 2)          
A_tP16_138_j-001                                 tP16      P4_{2}/ncm          A18       Cl2                 #138      Cl2 (part 1)        
A_tI2_139_a-001                                  tI2       I4/mmm              A6        In                  #139      In (part 1)         
A_tI2_139_a-002                                  tI2       I4/mmm              A_a       Pa                  #139      alpha-Pa (part 1)    
A_tI4_139_e-001                                  tI4       I4/mmm              -         Si                  #139      Hypothetical BCT5 Si (part 1)
A_tI8_139_h-001                                  tI8       I4/mmm              -         C                   #139      Hypothetical Tetrahedrally Bonded Carbon with 4-Member Rings (part 1)
A_tI4_141_a-001                                  tI4       I4_{1}/amd          A5        Sn                  #141      betaSn (part 1)     
A_tI16_142_f-001                                 tI16      I4_{1}/acd          -         S                   #142      S-III (part 2)      
A_hP3_152_a-001                                  hP3       P3_{1}21            A8        Se                  #152      gammaSe (part 1)    
A_hP9_154_bc-001                                 hP9       P3_{2}21            -         S                   #154      S-II (part 2)       
A_hR105_166_bc9h4i-001                           hR105     R-3m                -         B                   #166      beta-Boron (part 1)  
A_hR12_166_2h-001                                hR12      R-3m                -         B                   #166      alpha-Boron (part 1) 
A_hR1_166_a-001                                  hR1       R-3m                A_i       Po                  #166      beta-Po (part 1)     
A_hR1_166_a-002                                  hR1       R-3m                A10       Hg                  #166      alpha-Hg (part 1)    
A_hR2_166_c-001                                  hR2       R-3m                A7        As                  #166      alpha-As (part 1)    
A_hR2_166_c-002                                  hR2       R-3m                -         C                   #166      Rhombohedral Graphite (part 1)
A_hR2_166_c-003                                  hR2       R-3m                -         O2                  #166      beta-O2 (part 1)     
A_hR3_166_ac-001                                 hR3       R-3m                C19       Sm                  #166      alpha-Sm (part 1)    
A_hP6_178_a-001                                  hP6       P6_{1}22            -         Sc                  #178      Sc-V (part 2)       
A_hP4_186_ab-001                                 hP4       P6_{3}mc            -         C                   #186      Buckled Graphite (part 1)
A_hP1_191_a-001                                  hP1       P6/mmm              A_f       gamma-HgSn6-10      #191      Simple Hexagonal Lattice (part 1)
A_hP2_194_c-001                                  hP2       P6_{3}/mmc          A3        Mg                  #194      Hexagonal Close Packed (part 1)
A_hP4_194_ac-001                                 hP4       P6_{3}/mmc          A3'       La                  #194      alpha-La (part 1)    
A_hP4_194_bc-001                                 hP4       P6_{3}/mmc          A9        C                   #194      Hexagonal Graphite structure (part 1)
A_hP4_194_f-001                                  hP4       P6_{3}/mmc          -         C                   #194      Lonsdaleite (part 1)
A_hP6_194_h-001                                  hP6       P6_{3}/mmc          -         C                   #194      Hypothetical Tetrahedrally Bonded Carbon with 3-Member Rings (part 1)
A_cP8_198_2a-001                                 cP8       P2_{1}3             -         N2                  #198      alpha-N2 (part 1)   
A_cF240_202_h2i-001                              cF240     Fm-3                -         C                   #202      FCC C60 Buckminsterfullerine (part 2)
A_cP240_205_10d-001                              cP240     Pa-3                -         C                   #205      Simple Cubic C60 Buckminsterfullerine (part 2)
A_cP8_205_c-001                                  cP8       Pa-3                -         N2                  #205      alpha-N2 (part 1)   
A_cI16_206_c-001                                 cI16      Ia-3                -         Si                  #206      BC8 (part 1)        
A_cP20_213_cd-001                                cP20      P4_{1}32            A13       Mn                  #213      beta-Mn (part 1)     
A_cI58_217_ac2g-001                              cI58      I-43m               A12       Mn                  #217      alpha-Mn (part 1)    
A_cI16_220_c-001                                 cI16      I-43d               -         Li                  #220      High Pressure cI16 Li (part 1)
A_cP1_221_a                                      cP1       Pm-3m               A_h       Po                  #221      alpha-Po (part 1)    
A_cP46_223_dik-001                               cP46      Pm-3n               -         Si                  #223      Si46 Clathrate (part 1)
A_cF4_225_a                                      cF4       Fm-3m               A1        Cu                  #225      Face-Centered Cubic (part 1), equivalent to 1,2,1.cub,2.cub [HTQC]
A_cF136_227_aeg-001                              cF136     Fd-3m               -         Si                  #227      Si34 Clathrate (part 1)
A_cF8_227_a                                      cF8       Fd-3m               A4        C                   #227      Diamond (part 1)    
A_cI2_229_a                                      cI2       Im-3m               A2        W                   #229      Body-Centered Cubic (part 1)

*****************************************************************************************************************************************************
ANRL - BINARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes (part)
----------------------------------------------------------------------------------------------------------------------------------------------------
AB2_aP12_1_4a_8a-001.{AB}                        aP12      P1                  -         FeS2                #1        FeS2 (part 1)       
A2B_aP6_2_2i_i-001.{AB}                          aP6       P-1                 -         P2I4                #2        P2I4 (part 1)       
A2B_aP6_2_aei_i-001.{AB}                         aP6       P-1                 -         H2S                 #2        H2S (part 2)        
A2B_mP12_3_bc3e_2e-001.{AB}                      mP12      P2                  -         SiO2                #3        SiO2 (part 1)       
A8B5_mP13_6_a7b_3a2b-001.{AB}                    mP13      Pm                  -         Mo8P5               #6        Mo8P5 (part 2)      
AB_mP4_6_2b_2a-001.AB                            mP4       Pm                  -         FeNi                #6        FeNi (part 2)       
A2B_mP12_7_4a_2a-001.{AB}                        mP12      Pc                  -         H2S                 #7        H2S IV (part 2)     
A2B_mP18_7_6a_3a-001.{AB}                        mP18      Pc                  -         As2Ba               #7        As2Ba (part 2)      
A3B_mP16_7_6a_2a-001.{AB}                        mP16      Pc                  -         epsilon-WO3         #7        epsilon-WO3 (part 2)
A9B2_mP22_7_9a_2a-001.{AB}                       mP22      Pc                  -         Rh2Ga9              #7        Rh2Ga9 (part 2)     
A2B_mC144_9_24a_12a-001.{AB}                     mC144     Cc                  -         SiO2                #9        Monoclinic Low Tridymite (part 1)
A5B3_mC32_9_5a_3a-001.{AB}                       mC32      Cc                  -         alpha-P3N5          #9        alpha-P3N5 (part 2) 
AB3_mC16_9_a_3a-001.{AB}                         mC16      Cc                  -         H3Cl                #9        H3Cl (part 2)       
A2B_mP6_10_mn_bg-001.{AB}                        mP6       P2/m                -         delta-PdCl2         #10       delta-PdCl2 (part 2)
AB3_mP16_10_mn_3m3n-001.{AB}                     mP16      P2/m                -         H3Cl                #10       H3Cl (part 2)       
AB_mP6_10_en_am-001.AB                           mP6       P2/m                -         LiSn                #10       LiSn (part 2)       
AB_mP4_11_e_e-001.{AB}                           mP4       P2_{1}/m            -         NiTi                #11       NiTi2 (part 1)      
A5B2_mC14_12_a2i_i-001.{AB}                      mC14      C2/m                -         Au5Mn2              #12       Au5Mn2 (part 1)     
AB2_mC6_12_a_i-001.{AB}                          mC6       C2/m                C34       AuTe2               #12       Calaverite (part 1) 
AB3_mC16_12_g_ij-001.{AB}                        mC16      C2/m                D0_15     AlCl3               #12       AlCl3 (part 1)      
A2B_mP12_13_2g_ef-001.{AB}                       mP12      P2/c                -         H2S                 #13       H2S (part 2)        
A2B_mP12_14_2e_e-001.{AB}                        mP12      P2_{1}/c            C43       ZrO2                #14       Baddeleyite (part 1)
A2B_mP12_14_2e_e-002.{AB}                        mP12      P2_{1}/c            -         O2W                 #14       O2W binary oxide (R. Friedrich), ICSD #8217
A2B_mP12_14_2e_e-003.{AB}                        mP12      P2_{1}/c            -         O2V                 #14       O2V binary oxide (R. Friedrich), ICSD #647610
A2B_mP6_14_e_a-001.{AB}                          mP6       P2_{1}/c            -         gamma-PdCl2         #14       gamma-PdCl2 (part 2)
A7B8_mP120_14_14e_16e-001.{AB}                   mP120     P2_{1}/c            -         alpha-C7H8          #14       alpha-Toluene (part 2)
A2B5_mC28_15_f_e2f-001.{AB}                      mC28      C2/c                -         B2Pd5               #15       B2Pd5 (part 1)      
A2B5_mC28_15_f_e2f-002.{AB}                      mC28      C2/c                -         O5Sb2               #15       O5Sb2 binary oxide (R. Friedrich), ICSD #1422      
A2B_mC48_15_ae3f_2f-001.{AB}                     mC48      C2/c                -         SiO2                #15       Coesite (part 1)    
AB3_mC16_15_e_cf-001.{AB}                        mC16      C2/c                -         H3Cl                #15       H3Cl (part 2)       
AB_mC8_15_c_e-001.{AB}                           mC8       C2/c                B26       CuO                 #15       Tenorite (part 1)   
A2B_oP12_17_abe_e-001.{AB}                       oP12      P222_{1}            -         alpha-Ag2Se         #17       alpha-Naumannite (part 2)
AB3_oP16_18_ab_3c-001.{AB}                       oP16      P2_{1}2_{1}2        D0_{17}   BaS3                #18       BaS3 (part 1)
A2B_oP12_19_2a_a-001.{AB}                        oP12      P2_{1}2_{1}2_{1}    -         Ag2Se               #19       Naumannite (part 1) 
AB3_oP16_19_a_3a-001.{AB}                        oP16      P2_{1}2_{1}2_{1}    -         H3Cl                #19       H3Cl (part 2)       
A2B_oC24_20_abc_c-001.{AB}                       oC24      C222_{1}            -         SiO2                #20       Orthorhombic Tridymite (part 1)
AB2_oC6_21_a_k-001.{AB}                          oC6       C222                -         Ta2H                #21       Ta2H (part 2)       
AB_oF8_22_a_c-001.AB                             oF8       F222                -         FeS                 #22       FeS (part 2)        
A3B_oI32_23_ij2k_k-001.{AB}                      oI32      I222                -         H3S                 #23       H3S (part 2)        
AB_oP2_25_b_a-001.AB                             oP2       Pmm2                -         CdTe                #25       High-Pressure CdTe (part 1)
A2B_oP12_26_abc_ab-001.{AB}                      oP12      Pmc2_{1}            -         H2S                 #26       H2S (part 2)        
A2B_oP12_26_abc_ab-002.{AB}                      oP12      Pmc2_{1}            -         beta-SeO2           #26       beta-SeO2 (part 2)  
A5B_oP24_26_3a3b2c_ab-001.{AB}                   oP24      Pmc2_{1}            -         TlP5                #26       TlP5 (part 2)       
AB2_oP24_28_acd_2c3d-001.{AB}                    oP24      Pma2                C46       AuTe2               #28       Krennerite (part 1) 
A2B_oP12_29_2a_a-001.{AB}                        oP12      Pca2_{1}            -         ZrO2                #29       ZrO2 (part 2)       
AB2_oP12_29_a_2a-001.{AB}                        oP12      Pca2_{1}            -         FeS2                #29       Pyrite (part 2)     
A2B3_oP40_33_4a_6a-001.{AB}                      oP40      Pna2_{1}            -         kappa-Al2O3         #33       kappa-alumina (part 2)
AB_oP8_33_a_a-001.{AB}                           oP8       Pna2_{1}            -         CoAs                #33       Modderite (part 1)  
AB2_oP6_34_a_c-001.{AB}                          oP6       Pnn2                -         FeSb2               #34       FeSb2 (part 2)      
A2B_oC12_36_2a_a-001.{AB}                        oC12      Cmc2_{1}            C24       HgBr2               #36       HgBr2 (part 1)      
AB_oC8_36_a_a-001.{AB}                           oC8       Cmc2_{1}            -         HCl                 #36       HCl (part 2)        
A2B_oC12_38_de_ab-001.{AB}                       oC12      Amm2                -         Au2V                #38       Au2V (part 1)       
A2B3_oC40_39_2d_2c2d-001.{AB}                    oC40      Aem2                -         Ta3S2               #39       Ta3S2 (part 2)      
AB3_oC16_40_b_3b-001.{AB}                        oC16      Ama2                -         CeTe3               #40       CeTe3 (part 2)      
AB2_oC24_41_2a_2b-001.{AB}                       oC24      Aea2                C_e       PdSn2               #41       PdSn2 (part 1)      
AB4_oC20_41_a_2b-001.{AB}                        oC20      Aea2                D1_c      PtSn4               #41       PtSn4 (part 1)      
A10B3_oF52_42_2abce_ab-001.{AB}                  oF52      Fmm2                -         W3O10               #42       W3O10 (part 2)      
AB_oF8_42_a_a-001.AB                             oF8       Fmm2                -         BN                  #42       BN (part 2)         
AB2_oF72_43_ab_3b-001.{AB}                       oF72      Fdd2                C44       GeS2                #43       GeS2 (part 1)       
AB_oI4_44_a_b-001.AB                             oI4       Imm2                -         AsGa                #44       High-pressure GaAs (part 1)
A5B2_oP14_49_dehq_ab-001.{AB}                    oP14      Pccm                -         beta-Ta2O5          #49       beta-Ta2O5 (part 2) 
AB_oP4_51_e_f-001.AB                             oP4       Pmma                B19       AuCd                #51       beta'-AuCd (part 1) 
A2B_oP24_52_2e_cd-001.{AB}                       oP24      Pnna                -         GaCl2               #52       GaCl2 (part 2)      
A3B2_oP20_52_de_cd-001.{AB}                      oP20      Pnna                -         Sr2Bi3              #52       Sr2Bi3 (part 2)     
A2B_oP24_55_2g2h_gh-001.{AB}                     oP24      Pbam                -         GeAs2               #55       GeAs2 (part 2)      
A3B5_oP16_55_ch_agh-001.{AB}                     oP16      Pbam                -         Rh5Ge3              #55       Rh5Ge3 (part 2)     
A3B2_oP20_56_ce_e-001.{AB}                       oP20      Pccn                D5_11     Sb2O3               #56       Sb2O3 (part 1)      
AB_oP8_57_d_d-001.{AB}                           oP8       Pbcm                -         TlF                 #57       TlF-II, 57 (part 1) 
A2B_oP6_58_g_a-001.{AB}                          oP6       Pnnm                C50       alpha-PdCl2         #58       alpha-PdCl2 (part 2)
AB2_oP6_58_a_g-001.{AB}                          oP6       Pnnm                C35       CaCl2               #58       Hydrophilite (part 1)
AB2_oP6_58_a_g-002.{AB}                          oP6       Pnnm                -         etaFe2C             #58       eta-Fe2C (part 1)   
AB2_oP6_58_a_g-003.{AB}                          oP6       Pnnm                C18       FeS2                #58       Marcasite (part 1)  
A3B_oP8_59_bf_a-001.{AB}                         oP8       Pmmn                D0_a      TiCu3               #59       betaTiCu3 (part 1)  
AB_oP4_59_a_b-001.{AB}                           oP4       Pmmn                -         CuTe                #59       Vulcanite (part 1)  
A2B3_oP20_60_d_cd-001.{AB}                       oP20      Pbcn                -         Rh2S3               #60       Rh2S3 (part 2)      
A3B_oP32_60_3d_d-001.{AB}                        oP32      Pbcn                -         WO3                 #60       WO3 (part 2)        
A7B8_oP120_60_7d_8d-001.{AB}                     oP120     Pbcn                -         beta-C7H8           #60       beta-Toluene (part 2)
A2B_oP24_61_2c_c-001.{AB}                        oP24      Pbca                C21       TiO2                #61       Brookite (part 1)   
AB_oP16_61_c_c-001.{AB}                          oP16      Pbca                B_e       CdSb                #61       CdSb (part 1)       
AB_oP48_61_3c_3c-001.{AB}                        oP48      Pbca                -         Benzene             #61       Benzene (part 2)    
A2B3_oP20_62_2c_3c-001.{AB}                      oP20      Pnma                D5_{10}   Cr3C2               #62       Tongbaite (part 2)  
A2B_oP12_62_2c_c-001.{AB}                        oP12      Pnma                C37       Co2Si               #62       Co2Si (part 1)      
A2B_oP12_62_2c_c-002.{AB}                        oP12      Pnma                C25       HgCl2               #62       HgCl2 (part 1)      
A2B_oP12_62_2c_c-003.{AB}                        oP12      Pnma                C23       PbCl2               #62       Cotunnite (part 1)  
A2B_oP12_62_2c_c-004.{AB}                        oP12      Pnma                C29       SrH2                #62       SrH2 (part 2)      
A3B2_oP20_62_3c_2c-001.{AB}                      oP20      Pnma                D5_8      Sb2S3               #62       Stibnite (part 1)   
A3B7_oP40_62_cd_3c2d-001.{AB}                    oP40      Pnma                D10_1     C3Cr7               #62       C3Cr7 (part 1)      
A3B_oP16_62_cd_c-001.{AB}                        oP16      Pnma                D0_{20}   epsilon-NiAl3       #62       epsilon-NiAl3 (part 2)
A4B_oP20_62_2cd_c-001.{AB}                       oP20      Pnma                -         MgB4                #62       MgB4 (part 1)       
AB3_oP16_62_c_3c-001.{AB}                        oP16      Pnma                D0_{8}    MoO3                #62       Molybdite (part 2)  
AB3_oP16_62_c_cd-001.{AB}                        oP16      Pnma                D0_11     Fe3C                #62       Cementite (part 1)  
AB_oP8_62_c_c-001.{AB}                           oP8       Pnma                B16       GeS                 #62       GeS (part 1)        
AB_oP8_62_c_c-002.{AB}                           oP8       Pnma                B31       MnP                 #62       MnP (part 1)        
AB_oP8_62_c_c-003.{AB}                           oP8       Pnma                B27       FeB                 #62       FeB (part 1)        
AB_oP8_62_c_c-004.{AB}                           oP8       Pnma                B29       SnS                 #62       SnS (part 1)        
AB_oP8_62_c_c-005.{AB}                           oP8       Pnma                B14       FeAs                #62       Westerveldite (part 2)
AB_oP8_62_c_c-006.{AB}                           oP8       Pnma                -         HgO                 #62       binary oxide (R. Friedrich), ICSD #14124 
A2B_oC12_63_2c_c-001.{AB}                        oC12      Cmcm                C49       ZrSi2               #63       ZrSi2 (part 1)      
A6B_oC28_63_efg_c-001.{AB}                       oC28      Cmcm                D2_{h}    MnAl6               #63       MnAl6 (part 2)      
AB_oC8_63_c_c-001.{AB}                           oC8       Cmcm                B33       CrB                 #63       CrB (part 1)        
AB_oC8_63_c_c-002.AB                             oC8       Cmcm                -         AB                  #63       sc derivative (Y. Lederer)
A2B_oC24_64_2f_f-001.{AB}                        oC24      Cmce                -         H2S                 #64       H2S (part 2)        
A3B5_oC16_65_ah_bej-001.{AB}                     oC16      Cmmm                -         Ga3Pt5              #65       Ga3Pt5 (part 1)     
AB3_oC8_65_a_bf-001.{AB}                         oC8       Cmmm                L1_3      CdPt3               #65       Predicted CdPt3 (part 1), fcc derivative (Y. Lederer), equivalent to 11,12 [HTQC]
AB_oC8_65_j_g-001.{AB}                           oC8       Cmmm                -         alpha-IrV           #65       alpha-IrV (part 1)   
A3B_oC64_66_gi2lm_2l-001.{AB}                    oC64      Cccm                -         H3S                 #66       H3S (part 2)        
A3B_oC64_66_kl2m_bdl-001.{AB}                    oC64      Cccm                -         beta-ThI3           #66       beta-ThI3 (part 2)  
AB_oC8_67_a_g-001.{AB}                           oC8       Cmme                -         alpha-FeSe          #67       alpha-FeSe (part 2) 
AB_oC8_67_a_g-002.{AB}                           oC8       Cmme                -         alpha-PbO           #67       alpha-PbO (part 2)  
AB4_oC20_68_a_i-001.{AB}                         oC20      Ccce                -         PdSn4               #68       PdSn4 (part 2)      
AB_oF8_69_a_b-001.AB                             oF8       Fmmm                B24       TlF                 #69       TlF (part 1)        
A2B_oF24_70_e_a-001.{AB}                         oF24      Fddd                C54       TiSi2               #70       TiSi2 (part 1)      
AB2_oF48_70_f_fg-001.{AB}                        oF48      Fddd                D1_{f}    Mn2B                #70       Mn2B (part 2)       
A4B3_oI14_71_gh_cg-001.{AB}                      oI14      Immm                D7_{b}    Ta3B4               #71       Ta3B4 (part 2)      
AB2_oI6_71_a_g-001.{AB}                          oI6       Immm                -         MoPt2               #71       MoPt2 (part 1), fcc derivative (Y. Lederer), equivalent to 7,8 [HTQC]     
AB2_oI6_71_a_i-001.{AB}                          oI6       Immm                -         ReSi2               #71       ReSi2 (part 1)      
A2B_oI12_72_j_a-001.{AB}                         oI12      Ibam                C42       SiS2                #72       SiS2 (part 1)       
A2B_oI12_74_h_e-001.{AB}                         oI12      Imma                -         KHg2                #74       KHg2 (part 2)       
A4B_oI20_74_beh_e-001.{AB}                       oI20      Imma                D1_{b}    Al4U                #74       Al4U (part 2)       
A3B7_tP40_76_3a_7a-001.{AB}                      tP40      P4_{1}              -         Cs3P7               #76       Cs3P7 (part 2)      
A2B_tP48_77_8d_4d-001.{AB}                       tP48      P4_{2}              -         H2S                 #77       H2S III (part 2)    
AB2_tI48_80_2b_4b-001.{AB}                       tI48      I4_{1}              -         beta-NbO2           #80       beta-NbO2 (part 2)  
AB2_tP12_81_adg_2h-001.{AB}                      tP12      P-4                 -         GeSe2               #81       GeSe2 (part 2)      
A3B_tI32_82_3g_g-001.{AB}                        tI32      I-4                 D0_{e}    Ni3P                #82       Ni3P (part 2)       
A3B2_tP10_83_adk_j-001.{AB}                      tP10      P4/m                -         Ti2Ge3              #83       Ti2Ge3 (part 2)     
AB_tP16_84_cej_k-001.{AB}                        tP16      P4_{2}/m            B34       PdS                 #84       PdS (part 1)        
A2B_tP30_85_ab2g_cg-001.{AB}                     tP30      P4/n                -         SrBr2               #85       SrBr2 (part 2)      
AB3_tP32_86_g_3g-001.{AB}                        tP32      P4_{2}/n            -         Ti3P                #86       Ti3P (part 2)       
A4B5_tI18_87_h_ah-001.{AB}                       tI18      I4/m                -         Ti5Te4              #87       Ti5Te4 (part 1)     
AB4_tI10_87_a_h-001.{AB}                         tI10      I4/m                D1_a      Ni4Mo               #87       Ni4Mo (part 1)      
A4B_tI20_88_f_a-001.{AB}                         tI20      I4_{1}/a            -         ThCl4               #88       ThCl4 (part 2)      
AB2_tI96_88_2f_4f-001.{AB}                       tI96      I4_{1}/a            -         alpha-NbO2          #88       alpha-NbO2 (part 2) 
A2B_tP12_92_b_a-001.{AB}                         tP12      P4_{1}2_{1}2        -         SiO2                #92       alpha-Cristobalite (part 1)
A2B_tP12_92_b_a-002.{AB}                         tP12      P4_{1}2_{1}2        -         O2Ti                #92       O2Ti binary oxide (R. Friedrich), ICSD #161691
A2B_tP36_96_3b_ab-001.{AB}                       tP36      P4_{3}2_{1}2        -         SiO2                #96       Keatite (part 1)    
A2B_tI12_98_f_a-001.{AB}                         tI12      I4_{1}22            -         CdAs2               #98       CdAs2 (part 2)      
A2B3_tP20_102_2c_b2c-001.{AB}                    tP20      P4_{2}nm            -         Gd3Al2              #102      Gd3Al2 (part 2)     
AB4_tP10_103_a_d-001.{AB}                        tP10      P4cc                -         NbTe4               #103      NbTe4 (part 2)      
A5B7_tI24_107_ac_abd-001.{AB}                    tI24      I4mm                -         Co5Ge7              #107      Co5Ge7 (part 2)     
AB_tI4_107_a_a-001.AB                            tI4       I4mm                -         GeP                 #107      GeP (part 2)        
A3B5_tI32_108_ac_a2c-001.{AB}                    tI32      I4cm                -         Sr5Si3              #108      Sr5Si3 (part 2)     
AB_tI8_109_a_a-001.AB                            tI8       I4_{1}md            -         NbAs                #109      NbAs (part 2)       
A2B_tP12_111_2n_adf-001.{AB}                     tP12      P-42m               -         MnF2                #111      MnF2 (part 2)       
AB_tP8_111_n_n-001.AB                            tP8       P-42m               -         VN                  #111      VN (part 2)
AB3_tP8_113_a_ce-001.{AB}                        tP8       P-42_{1}m           -         BaS3                #113      BaS3 (part 1)
A3B_tP32_114_3e_e-001.{AB}                       tP32      P-42_{1}c           -         SeO3                #114      SeO3 (part 2)       
A4B_tP10_114_e_a-001.{AB}                        tP10      P-42_{1}c           -         Pd4Se               #114      Pd4Se (part 2)      
A2B3_tP5_115_g_ag-001.{AB}                       tP5       P-4m2               -         Rh3P2               #115      Rh3P2 (part 2)      
AB2_tP12_115_j_egi-001.{AB}                      tP12      P-4m2               -         HgI2                #115      HgI2 (part 2)       
A2B3_tP20_116_bci_fj-001.{AB}                    tP20      P-4c2               -         Ru2Sn3              #116      Ru2Sn3 (part 2)     
A2B3_tP20_117_i_adgh-001.{AB}                    tP20      P-4b2               D5_{12}   beta-Bi2O3          #117      beta-Bi2O3 (part 2)
A3B_tP16_118_ei_f-001.{AB}                       tP16      P-4n2               -         RuIn3               #118      RuIn3 (part 2)      
A5B3_tP32_118_g2i_aceh-001.{AB}                  tP32      P-4n2               -         Ir3Ga5              #118      Ir3Ga5 (part 2)     
A3B_tI24_119_b2i_af-001.{AB}                     tI24      I-4m2               -         RbGa3               #119      RbGa3 (part 2)      
AB_tI4_119_c_a-001.AB                            tI4       I-4m2               -         GaSb                #119      GaSb (part 2)       
AB3_tP4_123_a_ce-001.{AB}                        tP4       P4/mmm              L6_0      CuTi3               #123      CuTi3 (part 1)      
AB_tP2_123_a_d-001.AB                            tP2       P4/mmm              L1_0      CuAu                #123      CuAu (part 1), fcc derivative (Y. Lederer), equivalent to 3 [HTQC] 
AB4_tP10_124_a_m-001.{AB}                        tP10      P4/mcc              -         NbTe4               #124      NbTe4 (part 2)      
A4B_tP10_125_m_a-001.{AB}                        tP10      P4/nbm              -         PtPb4               #125      PtPb4 (part 2)      
A2B3_tP10_127_g_ah-001.{AB}                      tP10      P4/mbm              D5_a      Si2U3               #127      Si2U3 (part 1)      
A4B_tP20_127_ehj_g-001.{AB}                      tP20      P4/mbm              D1_{e}    ThB4                #127      ThB4 (part 2)       
A2B_tP6_129_ac_c-001.{AB}                        tP6       P4/nmm              C38       Cu2Sb               #129      Cu2Sb (part 1)      
AB_tP4_129_a_c-001.{AB}                          tP4       P4/nmm              B10       PbO                 #129      PbO (part 1)        
AB_tP4_129_c_c-001.{AB}                          tP4       P4/nmm              B11       gamma-CuTi          #129      gamma-CuTi (part 1)  
AB_tP4_129_c_c-002.AB                            tP4       P4/nmm              -         AB                  #129      fcc derivative (Y. Lederer), equivalent to 14 [HTQC]
A5B3_tP32_130_cg_cf-001.{AB}                     tP32      P4/ncc              -         Ba5Si3              #130      Ba5Si3 (part 2)     
AB_tP4_131_c_e-001.{AB}                          tP4       P4_{2}/mmc          B17       PtS                 #131      PtS (part 1)        
AB_tP4_131_c_e-002.{AB}                          tP4       P4_{2}/mmc          -         PdO                 #131      PdO binary oxide (R. Friedrich), ICSD #26598        
AB3_tP32_133_h_i2j-001.{AB}                      tP32      P4_{2}/nbc          -         beta-V3S            #133      beta-V3S (part 2)   
A2B_tP24_135_gh_h-001.{AB}                       tP24      P4_{2}/mbc          C47       SeO2                #135      Downeyite (part 2)  
A2B_tP6_136_f_a-001.{AB}                         tP6       P4_{2}/mnm          C4        TiO2                #136      CrSi2 (part 1)      
AB_tP8_136_g_f-001.{AB}                          tP8       P4_{2}/mnm          -         beta-BeO            #136      betaBeO (part 1)    
A2B3_tP40_137_cdf_3g-001.{AB}                    tP40      P4_{2}/nmc          D5_{9}    Zn3P2               #137      Zn3P2 (part 2)      
A2B_tP6_137_d_a-001.{AB}                         tP6       P4_{2}/nmc          -         ZrO2                #137      ZrO2 (part 2)       
AB2_tP6_137_a_d-001.{AB}                         tP6       P4_{2}/nmc          C13       HgI2                #137      HgI2 (part 2)       
A12B_tI26_139_fij_a-001.{AB}                     tI26      I4/mmm              D2_b      Mn12Th              #139      Mn12Th (part 1)     
A2B_tI6_139_d_a-001.{AB}                         tI6       I4/mmm              L\'2      ThH2                #139      ThH2 (part 1)       
A3B_tI16_139_cde_e-001.{AB}                      tI16      I4/mmm              D0_23     Al3Zr               #139      Al3Zr (part 1)      
A3B_tI8_139_bd_a-001.{AB}                        tI8       I4/mmm              D0_22     Al3Ti               #139      Al3Ti (part 1)      
A3B_tI8_139_bd_a-002.{AB}                        tI8       I4/mmm              -         AB3                 #139      fcc derivative (Y. Lederer) equivalent to 22,24 [HTQC]
A4B5_tI18_139_i_ah-001.{AB}                      tI18      I4/mmm              -         V4Zn5               #139      V4Zn5 (part 1)      
A4B_tI10_139_de_a-001.{AB}                       tI10      I4/mmm              D1_3      Al4Ba               #139      Al4Ba (part 1)      
A8B_tI18_139_hi_a-001.{AB}                       tI18      I4/mmm              -         Pt8Ti               #139      Pt8Ti (part 1)      
AB2_tI6_139_a_e-001.{AB}                         tI6       I4/mmm              C11_b     MoSi2               #139      MoSi2 (part 1)      
AB2_tI6_139_a_e-002.{AB}                         tI6       I4/mmm              -         AB2                 #139      fcc derivative (Y. Lederer), equivalent to 5,6,6b [HTQC]
AB2_tI6_139_a_e-003.{AB}                         tI6       I4/mmm              -         BaO2                #139      BaO2 binary oxide (R. Friedrich), ICSD #24248 
AB2_tI6_139_a_e-004.{AB}                         tI6       I4/mmm              -         KO2                 #139      KO2 binary oxide (R. Friedrich), ICSD #38245, equivalent to 250,251,1202.AB/BA,1211.AB/BA [HTQC]
AB_tI8_139_e_e-001.{AB}                          tI8       I4/mmm              D3_{1}    Hg2Cl2              #139      Calomel (part 2)    
AB_tI8_139_e_e-002.AB                            tI8       I4/mmm              -         AB                  #139      sc derivative (Y. Lederer)
A2B_tI12_140_h_a-001.{AB}                        tI12      I4/mcm              C16       Al2Cu               #140      Khatyrkite (part 1) 
A3B5_tI32_140_ah_bk-001.{AB}                     tI32      I4/mcm              D8_{m}    W5Si3               #140      W5Si3 (part 2)      
A3B5_tI32_140_ah_cl-001.{AB}                     tI32      I4/mcm              D8_{l}    Cr5B3               #140      Cr5B3 (part 2)      
AB3_tI16_140_b_ah-001.{AB}                       tI16      I4/mcm              D0_c      SiU3                #140      SiU3 (part 1)       
AB_tI16_140_ab_h-001.{AB}                        tI16      I4/mcm              B37       SeTl                #140      SeTl (part 1)       
A2B3_tI80_141_ceh_3h-001.{AB}                    tI80      I4_{1}/amd          -         In2S3               #141      betaIn2S3 (part 1)  
A2B_tI12_141_e_a-001.{AB}                        tI12      I4_{1}/amd          C5        TiO2                #141      Anatase (part 1)    
A2B_tI12_141_e_a-002.{AB}                        tI12      I4_{1}/amd          C_{c}     alpha-ThSi2         #141      alpha-ThSi2 (part 2)
A2B_tI24_141_2e_e-001.{AB}                       tI28      I4_{1}/amd          -         Ga2Hf               #141      Ga2Hf (part 1)      
A3B4_tI28_141_ad_h-001.{AB}                      tI28      I4_{1}/amd          -         Mn3O4               #141      Hausmannite (part 1)
AB_tI16_141_e_e-001.{AB}                         tI16      I4_{1}/amd          B_g       MoB                 #141      MoB (part 1)        
AB_tI8_141_a_b-001.AB                            tI8       I4_{1}/amd          "40"      NbP                 #141      NbP (part 1)        
AB_tI8_141_a_b-002.AB                            tI8       I4_{1}/amd          -         AB                  #141      fcc derivative (Y. Lederer), equivalent to 23 [HTQC]
AB2_hP12_143_cd_ab2d-001.{AB}                    hP12      P3                  -         MoS2                #143      MoS2 (part 2)       
A4B_hP15_144_4a_a-001.{AB}                       hP15      P3_{1}              -         IrGe4               #144      IrGe4 (part 2)      
AB_hP6_144_a_a-001.AB                            hP6       P3_{1}              -         ZnTe                #144      ZnTe (part 2)
A2B_hP9_147_g_ad-001.{AB}                        hP9       P-3                 B_b       AgZn                #147      zetaAgZn (part 1)   
A2B_hR18_148_2f_f-001.{AB}                       hR18      R-3                 -         beta-PdCl2          #148      beta-PdCl2 (part 2) 
AB3_hR8_148_c_f-001.{AB}                         hR8       R-3                 D0_5      BiI3                #148      BiI3 (part 1)       
AB_hR16_148_cf_cf-001.{AB}                       hR16      R-3                 -         C8H8                #148      Solid Cubane (part 1)
AB_hR26_148_b2f_a2f-001.AB                       hR26      R-3                 -         PdAl                #148      PdAl (part 1)       
AB3_hP24_149_acgi_3l-001.{AB}                    hP24      P312                -         Ti3O                #149      Ti3O (part 2)       
A2B_hP9_150_ef_bd-001.{AB}                       hP9       P321                C22       Fe2P                #150      Original Fe2P (part 1)
A3B_hP24_151_3c_2a-001.{AB}                      hP24      P3_{1}12            D0_4      CrCl3               #151      CrCl3 (part 1)      
A2B_hP9_152_c_a-001.{AB}                         hP9       P3_{1}21            -         SiO2                #152      alphaQuartz (part 1)
A3B_hP24_153_3c_2b-001.{AB}                      hP24      P3_{2}12            -         CrCl3               #153      CrCl3 (part 2)      
AB_hP6_154_a_b-001.{AB}                          hP6       P3_{2}21            B9        HgS                 #154      Cinnabar (part 1)   
A3B2_hR5_155_e_c-001.{AB}                        hR5       R32                 D5_e      Ni3S2               #155      Hazelwoodite (part 1)
AB3_hR8_155_c_de-001.{AB}                        hR8       R32                 D0_14     AlF3                #155      AlF3 (part 1)       
AB2_hP9_156_b2c_3a2bc-001.{AB}                   hP9       P3m1                -         CdI2                #156      CdI2 (part 2)       
AB_hP12_156_2ab3c_2ab3c-001.{AB}                 hP12      P3m1                -         CuI                 #156      CuI (part 2)        
AB_hP4_156_ab_ab-001.{AB}                        hP4       P3m1                -         beta-CuI            #156      beta-CuI (part 2)
A3B_hP8_158_d_a-001.{AB}                         hP8       P3c1                -         beta-RuCl3          #158      beta-RuCl3 (part 2) 
A2B3_hP20_159_bc_2c-001.{AB}                     hP20      P31c                -         Bi2O3               #159      Bi2O3 (part 2)      
A4B3_hP28_159_ab2c_2c-001.{AB}                   hP28      P31c                -         alpha-Si3N4         #159      Nierite (part 2)    
A3B_hR4_160_b_a-001.{AB}                         hR4       R3m                 -         H3S                 #160      H3S (part 2)        
A8B5_hR26_160_a3bc_a3b-001.{AB}                  hR26      R3m                 D8_{10}   Cr5Al8              #160      Cr5Al8 (part 2)
AB_hR10_160_5a_5a-001.AB                         hR10      R3m                 B7        SiC                 #160      Moissanite-15R (part 2)
AB_hR6_160_3a_3a-001.AB                          hR6       R3m                 -         CSi                 #160      Moissanite 9R structure (part 1)
AB_hR6_160_b_b-001.{AB}                          hR6       R3m                 B13       NiS                 #160      Millerite (part 1)  
AB2_hP9_162_ad_k-001.{AB}                        hP9       P-31m               -         NV2                 #162      betaV2N (part 1)    
A2B3_hP5_164_d_ad-001.{AB}                       hP5       P-3m1               D5_{2}    La2O3               #164      La2O3 (part 2)      
A3B2_hP5_164_ad_d-001.{AB}                       hP5       P-3m1               D5_13     Al3Ni2              #164      Al3Ni2 (part 1)     
AB2_hP3_164_a_d-001.{AB}                         hP3       P-3m1               C6        CdI2                #164      omega Phase (part 1)
AB2_hP9_164_bd_c2d-001.{AB}                      hP9       P-3m1               -         deltaH^II-NW2       #164      deltaH^II-NW2 (part 2)
A3B_hP24_165_adg_f-001.{AB}                      hP24      P-3c1               -         H3Ho                #165      H3Ho (part 1)       
A3B_hP24_165_bdg_f-001.{AB}                      hP24      P-3c1               D0_{21}   Cu3P                #165      Cu3P (part 2)       
A2B3_hR5_166_c_ac-001.{AB}                       hR5       R-3m                C33       Bi2Te3              #166      Bi2Te3 (part 1)     
A4B3_hR7_166_2c_ac-001.{AB}                      hR7       R-3m                D7_{1}    Al4C3               #166      Al4C3 (part 2)      
A5B2_hR7_166_a2c_c-001.{AB}                      hR7       R-3m                D8_i      B5Mo2               #166      Mo2B5 (part 1)      
A6B_hR7_166_g_a-001.{AB}                         hR7       R-3m                -         CaC6                #166      CaC6 (part 1)       
A7B6_hR13_166_ah_3c-001.{AB}                     hR13      R-3m                D8_5      Fe7W6               #166      Fe7W6 mu-phase (part 1)
AB_hR2_166_a_b-001.AB                            hR2       R-3m                L1_1      CuPt                #166      CuPt (part 1), fcc derivative (Y. Lederer), equivalent to 4 [HTQC] 
A2B3_hR10_167_c_e-001.{AB}                       hR10      R-3c                D5_1      Al2O3               #167      Corundum (part 1)   
A2B3_hR10_167_c_e-002.{AB}                       hR10      R-3c                -         O3V2                #167      O3V2 binary oxide (R. Friedrich), ICSD #94768
A2B3_hP30_169_2a_3a-001.{AB}                     hP30      P6_{1}              -         alpha-Al2S3         #169      alpha-Al2S3 (part 2)
A2B3_hP30_170_2a_3a-001.{AB}                     hP30      P6_{5}              -         Al2S3               #170      Al2S3 (part 2)      
A3B_hP8_173_c_b-001.{AB}                         hP8       P6_{3}              -         PI3                 #173      PI3 (part 2)        
A4B3_hP14_173_bc_c-001.{AB}                      hP14      P6_{3}              -         beta-Si3N4          #173      beta-Si3N4 (part 2) 
A3B2_hP10_176_h_bc-001.{AB}                      hP10      P6_{3}/m            -         Er3Ru2              #176      Er3Ru2 (part 2)
A3B_hP8_176_h_c-001.{AB}                         hP8       P6_{3}/m            -         UCl3                #176      UCl3 (part 2)
A2B_hP36_177_j2lm_n-001.{AB}                     hP36      P622                -         SiO2                #177      SiO2 (part 2)       
AB3_hP24_178_b_ac-001.{AB}                       hP24      P6_{1}22            -         AuF3                #178      AuF3 (part 2)       
AB3_hP24_179_b_ac-001.{AB}                       hP24      P6_{5}22            -         AuF3                #179      AuF3 (part 2)       
A2B_hP18_180_fi_bd-001.{AB}                      hP18      P6_{2}22            C_a       Mg2Ni               #180      Mg2Ni (part 1)      
A2B_hP9_180_j_c-001.{AB}                         hP9       P6_{2}22            C8        SiO2                #180      beta-Quartz (part 1)
AB2_hP9_180_d_j-001.{AB}                         hP9       P6_{2}22            C40       CrSi2               #180      CrSi2 (part 1)      
A2B_hP9_181_j_c-001.{AB}                         hP9       P6_{4}22            -         beta-SiO2           #181      beta-SiO2 (part 2)  
AB3_hP8_182_c_g-001.{AB}                         hP8       P6_{3}22            -         Fe3C                #182      Bainite (part 1)    
AB_hP6_183_c_ab-001.{AB}                         hP6       P6mm                -         CrFe3NiSn5          #183      CrFe3NiSn5 (part 2) 
A3B_hP24_185_ab2c_c-001.{AB}                     hP24      P6_{3}cm            -         Cu3P                #185      Cu3P (part 2)       
A3B_hP8_185_c_a-001.{AB}                         hP8       P6_{3}cm            -         beta-RuCl3          #185      beta-RuCl3 (part 2) 
AB3_hP24_185_c_ab2c-001.{AB}                     hP24      P6_{3}cm            -         Na3As               #185      Na3As (part 2)      
A3B7_hP20_186_c_b2c-001.{AB}                     hP20      P6_{3}mc            D10_{2}   Fe3Th7              #186      Fe3Th7 (part 2)     
AB_hP12_186_a2b_a2b-001.AB                       hP12      P6_{3}mc            B6        SiC                 #186      Moissanite-6H SiC (part 1)
AB_hP4_186_b_a-001.{AB}                          hP4       P6_{3}mc            B12       BN                  #186      Original BN (part 1)
AB_hP4_186_b_b-001.AB                            hP4       P6_{3}mc            B4        ZnS                 #186      Wurtzite (part 1)   
AB_hP8_186_ab_ab-001.AB                          hP8       P6_{3}mc            B5        SiC                 #186      Moissanite-4H SiC (part 1)
AB3_hP4_187_e_fh-001.{AB}                        hP4       P-6m2               -         Re3N                #187      Re3N (part 2)       
AB_hP2_187_d_a-001.AB                            hP2       P-6m2               B_h       WC                  #187      Tungsten Carbide (part 1)
A2B_hP9_189_fg_bc-001.{AB}                       hP9       P-62m               C22       Fe2P                #189      Revised Fe2P (part 1)
A2B_hP18_190_gh_bf-001.{AB}                      hP18      P-62c               -         Li2Sb               #190      Li2Sb (part 2)      
A5B3_hP16_190_bdh_g-001.{AB}                     hP16      P-62c               -         alpha-Sm3Ge5        #190      alpha-Sm3Ge5 (part 2)
AB_hP24_190_i_afh-001.{AB}                       hP24      P-62c               -         FeS                 #190      Troilite (part 2)   
A2B_hP6_191_h_e-001.{AB}                         hP6       P6/mmm              C_h       Cu2Te               #191      Cu2Te (part 1)      
A3B_hP4_191_bc_a-001.{AB}                        hP4       P6/mmm              -         Li3N                #191      Li3 Ni (part 1)     
AB2_hP3_191_a_d-001.{AB}                         hP3       P6/mmm              C32       AlB2                #191      Hexagonal omega (part 1)
AB5_hP6_191_a_cg-001.{AB}                        hP6       P6/mmm              D2_d      CaCu5               #191      CaCu5 (part 1)      
AB_hP6_191_f_ad-001.{AB}                         hP6       P6/mmm              B35       CoSn                #191      CoSn (part 1)       
AB2_hP72_192_m_j2kl-001.{AB}                     hP72      P6/mcc              -         AlPO4               #192      AlPO4 (part 2)      
A5B3_hP16_193_dg_g-001.{AB}                      hP16      P6_{3}/mcm          D8_{8}    Mn5Si3              #193      Mavlyanovite (part 2)
A2B_hP12_194_cg_f-001.{AB}                       hP12      P6_{3}/mmc          C10       SiO2                #194      beta-Tridymite (part 1)
A3B_hP16_194_gh_ac-001.{AB}                      hP16      P6_{3}/mmc          D0_{24}   Ni3Ti               #194      Ni3Ti (part 2)      
A3B_hP8_194_h_c-001.{AB}                         hP8       P6_{3}/mmc          D0_19     Ni3Sn               #194      Ni3Sn (part 1)      
A5B2_hP14_194_abdf_f-001.{AB}                    hP14      P6_{3}/mmc          D8_h      W2B5                #194      W2B5 (part 1)       
A5B2_hP28_194_ahk_ch-001.{AB}                    hP28      P6_{3}/mmc          D8_{11}   Co2Al5              #194      Co2Al5 (part 2)     
AB2_hP12_194_f_ah-001.{AB}                       hP12      P6_{3}/mmc          C14       MgZn2               #194      MgZn2 Hexagonal Laves (part 1)
AB2_hP24_194_ef_fgh-001.{AB}                     hP24      P6_{3}/mmc          C36       MgNi2               #194      MgNi2 Hexagonal Laves (part 1)
AB2_hP6_194_b_f-001.{AB}                         hP6       P6_{3}/mmc          -         CaIn2               #194      CaIn2 (part 1)      
AB2_hP6_194_c_ad-001.{AB}                        hP6       P6_{3}/mmc          B8_2      Ni2In               #194      Ni2In (part 1)      
AB2_hP6_194_c_f-001.{AB}                         hP6       P6_{3}/mmc          C7        MoS2                #194      Molybdenite (part 1)
AB3_hP8_194_c_bf-001.{AB}                        hP8       P6_{3}/mmc          D0_18     Na3As               #194      Na3As (part 1)      
AB_hP12_194_af_bf-001.{AB}                       hP12      P6_{3}/mmc          -         CMo                 #194      CMo structure (part 1)
AB_hP12_194_df_ce-001.{AB}                       hP12      P6_{3}/mmc          B18       CuS                 #194      Covellite (part 1)  
AB_hP4_194_c_a-001.{AB}                          hP4       P6_{3}/mmc          B8_1      NiAs                #194      NiAs (part 1)       
AB_hP4_194_c_d-001.AB                            hP4       P6_{3}/mmc          B_k       BN                  #194      BN (part 1)         
AB_hP8_194_ad_f-001.{AB}                         hP8       P6_{3}/mmc          B_i       AsTi                #194      AsTi (part 1)       
A4B_cI40_197_cde_c-001.{AB}                      cI40      I23                 -         Ga4Ni               #197      Ga4Ni3 (part 1)     
A3B_cP16_198_b_a-001.{AB}                        cP16      P2_{1}3             D1        NH3                 #198      Ammonia (part 1)    
AB_cP8_198_a_a-001.{AB}                          cP8       P2_{1}3             B21       CO                  #198      alphaCO (part 1)    
AB_cP8_198_a_a-002.AB                            cP8       P2_{1}3             B20       FeSi                #198      FeSi (part 1)       
AB_cI16_199_a_a-001.AB                           cI16      I2_{1}3             B_a       CoU                 #199      CoU (part 1)        
A2B11_cP39_200_f_aghij-001.{AB}                  cP39      Pm-3                D8_{c}    Mg2Zn11             #200      Mg2Zn11 (part 2)    
A12B_cI26_204_g_a-001.{AB}                       cI26      Im-3                -         Al12W               #204      Al12W (part 1)      
A3B_cI32_204_g_c-001.{AB}                        cI32      Im-3                D0_2      CoAs3               #204      Skutterudite (part 1)
AB2_cP12_205_a_c-001.{AB}                        cP12      Pa-3                C2        FeS2                #205      Pyrite (part 1)     
AB2_cP12_205_a_c-002.{AB}                        cP12      Pa-3                -         NaO2                #205      binary oxide (R. Friedrich), ICSD #87178, equivalent to 1214 [HTQC]
AB_cP16_205_c_c-001.{AB}                         cP16      Pa-3                -         CuCl                #205      SC16 (part 1)       
A17B15_cP64_207_acfk_eij-001.{AB}                cP64      P432                -         Pd17Se15            #207      Pd17Se15 (part 2)   
A3B_cP16_208_j_b-001.{AB}                        cP16      P4_{2}32            -         PH3                 #208      PH3 (part 2)        
A2B_cI72_211_hi_i-001.{AB}                       cI72      I432                -         SiO2                #211      SiO2 (part 2)       
A2B_cP12_212_c_a-001.{AB}                        cP12      P4_{3}32            -         SrSi2               #212      SrSi2 (part 2)      
A4B9_cP52_215_ei_3efgi-001.{AB}                  cP52      P-43m               D8_{3}    gamma-Cu9Al4        #215      gamma-brass (part 2)
AB4_cP5_215_a_e-001.{AB}                         cP5       P-43m               -         Fe4C                #215      Fe4C (part 1)       
AB5_cF24_216_a_ce-001.{AB}                       cF24      F-43m               C15_b     AuBe5               #216      AuBe5 (part 1)      
AB_cF8_216_c_a.AB                                cF8       F-43m               B3        ZnS                 #216      ZnS (part 1)       
A4B_cI10_217_c_a-001.{AB}                        cI10      I-43m               -         SiF4                #217      SiF4 (part 1)       
A5B8_cI52_217_ce_cg-001.{AB}                     cI52      I-43m               -         Cu5Zn8              #217      gammaBrass (part 1) 
A15B4_cI76_220_ae_c-001.{AB}                     cI76      I-43d               D8_{6}    Cu15Si4             #220      Cu15Si4 (part 2)    
A3B2_cI40_220_d_c-001.{AB}                       cI40      I-43d               D5_c      Pu2C3               #220      Pu2C3 (part 1)      
A4B3_cI28_220_c_a-001.{AB}                       cI28      I-43d               D7_{3}    Th3P4               #220      Th3P4 (part 2)      
A3B_cP4_221_d_a.{AB}                             cP4       Pm-3m               D0_9      ReO3                #221      alpha-ReO3 (part 1)  
A6B_cP7_221_f_a-001.{AB}                         cP7       Pm-3m               D2_1      CaB6                #221      CaB6 (part 1)       
AB11_cP36_221_c_agij-001.{AB}                    cP36      Pm-3m               D2_e      BaHg11              #221      BaHg11 (part 1)     
AB3_cP4_221_a_c.{AB}                             cP4       Pm-3m               L1_2      Cu3Au               #221      Cu3Au (part 1), fcc derivative (Y. Lederer), equivalent to 25,26 [HTQC] 
AB_cP2_221_b_a.AB                                cP2       Pm-3m               B2        CsCl                #221      CsCl (part 1), equivalent to 61,1026.AB,1205.AB [HTQC]
AB_cP6_221_c_d.AB                                cP6       Pm-3m               -         NbO                 #221      NbO (part 1)        
A3B_cP8_223_c_a.{AB}                             cP8       Pm-3n               A15       Cr3Si               #223      Cr3Si (part 1)      
A2B_cP6_224_b_a.{AB}                             cP6       Pn-3m               C3        Cu2O                #224      Cuprite (part 1)    
A12B_cF52_225_i_a-001.{AB}                       cF52      Fm-3m               D2_f      UB12                #225      UB12 (part 1)       
A23B6_cF116_225_bd2f_e-001.{AB}                  cF116     Fm-3m               D8_{a}    Th6Mn23             #225      Th6Mn23 (part 2)    
A6B23_cF116_225_e_acfh-001.{AB}                  cF116     Fm-3m               D8_4      Cr23C6              #225      Cr23C6 (part 1)     
A7B_cF32_225_bd_a.{AB}                           cF32      Fm-3m               -         Ca7Ge               #225      Ca7Ge (part 1)      
AB2_cF12_225_a_c.{AB}                            cF12      Fm-3m               C1        CaF2                #225      CaF2 (part 1)
AB3_cF16_225_a_bc.{AB}                           cF16      Fm-3m               D0_3      BiF3                #225      BiF3 (part 1)       
AB_cF8_225_a_b.AB                                cF8       Fm-3m               B1        NaCl                #225      Rock Salt (part 1), equivalent to 201,720.AB,720.BA,1009.AB,1200.AB [HTQC] 
AB13_cF112_226_a_bi-001.{AB}                     cF112     Fm-3c               D2_{3}    NaZn13              #226      NaZn13 (part 2)     
A2B_cF24_227_c_a.{AB}                            cF24      Fd-3m               C9        SiO2                #227      Ideal beta-Cristobalite (part 1)
A2B_cF24_227_d_a.{AB}                            cF24      Fd-3m               C15       Cu2Mg               #227      Cu2Mg Cubic Laves (part 1)
A3B4_cF56_227_ad_e-001.{AB}                      cF56      Fd-3m               D7_{2}    Co3O4               #227      Spinel (part 2)     
AB2_cF48_227_c_e-001.{AB}                        cF48      Fd-3m               -         CTi2                #227      CTi2 (part 1)       
AB2_cF96_227_e_cf-001.{AB}                       cF96      Fd-3m               -         NiTi2               #227      NiTi2 (part 1)      
AB_cF16_227_a_b.AB                               cF16      Fd-3m               B32       NaTl                #227      NaTl (part 1)       
A6B_cF224_228_h_c-001.{AB}                       cF224     Fd-3c               -         TeO6H6              #228      TeO6H6 (part 2)     
A2B7_cI54_229_e_afh-001.{AB}                     cI54      Im-3m               L2_2      Sb2Tl7              #229      L22 (part 1)        
A3B10_cI52_229_e_fh-001.{AB}                     cI52      Im-3m               D8_{1}    gamma-Fe3Zn10       #229      gamma-brass (part 2)
A3B_cI8_229_b_a.{AB}                             cI8       Im-3m               -         H3S                 #229      High Pressure H3S (part 1)
A4B3_cI14_229_c_b.{AB}                           cI14      Im-3m               -         Pt3O4               #229      Pt3O4 (part 1)      
A4B_cI10_229_c_a.{AB}                            cI10      Im-3m               -         beta-Hg4Pt          #229      beta-Hg4Pt (part 2) 
A7B3_cI40_229_df_e-001.{AB}                      cI40      Im-3m               D8_{f}    Ir3Ge7              #229      Ir3Ge7 (part 2)     
A4B3_cI112_230_af_g-001.{AB}                     cI112     Ia-3d               -         Ga4Ni3              #230      Ga4Ni3 (part 1)     

*****************************************************************************************************************************************************
ANRL - TERNARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes (part)
----------------------------------------------------------------------------------------------------------------------------------------------------
ABC2_aP16_1_4a_4a_8a-001.{ABC}                   aP16      P1                  -         AsKSe2              #1        AsKSe2 (part 1)     
A3BC_mC10_8_ab_a_a-001.{ABC}                     mC10      Cm                  -         Pb(Zr0.52Ti0.48)O3  #8        Monoclinic PZT [PbO3] (part 1)
ABC2_mP8_10_ac_eh_mn-001.ABC/BCA/CAB             mP8       P2/m                -         AuAgTe2             #10       Muthmannite (part 2)
ABC3_mP10_11_e_e_ef-001.{ABC}                    mP10      P2_{1}/m            G0_6      KClO3               #11       KClO3 (part 1)      
A7B2C2_mC22_12_aij_h_i-001.{ABC}                 mC22      C2/m                S2_{1}    [Sc,Y]2Si2O7        #12       Thortveitite (part 2)
ABC4_mP12_13_e_a_2g-001.{ABC}                    mP12      P2/c                E1_b      AgAuTe4             #13       Sylvanite (part 1)  
ABC4_oP12_16_ag_cd_2u-001.{ABC}                  oP12      P222                -         AlPS4               #16       AlPS4 (part 1)      
A2BC2_oF40_22_fi_ad_gh-001.{ABC}                 oF40      F222                -         CeRu2B2             #22       CeRu2B2 (part 2)    
ABC2_oI16_23_ab_i_k-001.{ABC}                    oI16      I222                -         NaFeS2              #23       NaFeS2 (part 2)     
ABC4_oI12_23_a_b_k-001.{ABC}                     oI12      I222                -         BPS4                #23       BPS4 (part 2)       
ABC_oP12_29_a_a_a-001.ABC/BCA/CAB                oP12      Pca2_{1}            -         CoAsS               #29       Cobaltite (part 2)  
A5B3C15_oP46_30_a2c_bc_a7c-001.{ABC}             oP46      Pnc2                -         Bi5Nb3O15           #30       Bi5Nb3O15 (part 2)  
ABC3_oP20_30_2a_c_3c-001.{ABC}                   oP20      Pnc2                -         CuBrSe3             #30       CuBrSe3 (part 2)    
AB3C4_oP16_31_a_ab_2ab-001.{ABC}                 oP16      Pmn2_{1}            H2_5      AsCu3S4             #31       Enargite (part 1)   
A13B2C2_oP34_32_a6c_c_c-001.{ABC}                oP34      Pba2                -         Re2O5[SO4]2         #32       Re2O52 (part 2)     
AB3C4_oP32_33_a_3a_4a-001.{ABC}                  oP32      Pna2_{1}            -         AsK3S4              #33       AsK3S4 (part 1)     
A2B8C_oP22_34_c_4c_a-001.{ABC}                   oP22      Pnn2                -         TiAl2Br8            #34       TiAl2Br8 (part 2)   
AB8C2_oC22_35_a_ab3e_e-001.{ABC}                 oC22      Cmm2                -         V2MoO8              #35       V2MoO8 (part 2)     
A2B5C2_oC36_37_d_c2d_d-001.{ABC}                 oC36      Ccc2                -         Li2Si2O5            #37       Li2Si2O5 (part 2)   
A2BC_oC8_38_e_a_b-001.{ABC}                      oC8       Amm2                -         C2CeNi              #38       C2CeNi (part 1)     
A9BC_oC44_39_3c3d_a_c-001.{ABC}                  oC44      Aem2                -         VPCl9               #39       VPCl9 (part 2)      
AB2C_oC16_40_a_2b_b-001.{ABC}                    oC16      Ama2                -         K2CdPb              #40       K2CdPb (part 2)     
A2BC2_oI20_45_c_b_c-001.ABC/BCA/CAB              oI20      Iba2                -         MnGa2Sb2            #45       MnGa2Sb2 (part 2)   
ABC_oI36_46_ac_bc_3b-001.{ABC}                   oI36      Ima2                -         TiFeSi              #46       TiFeSi (part 2)     
A2BC4_oP28_50_ij_ac_ijm-001.{ABC}                oP28      Pban                -         La2NiO4             #50       La2NiO4 (part 2)    
A3BC2_oP48_50_3m_m_2m-001.{ABC}                  oP48      Pban                -         alpha-Tl2TeO3       #50       alpha-Tl2TeO3 (part 2)
ABC2_oP16_53_h_e_gh-001.{ABC}                    oP16      Pmna                -         TaNiTe2             #53       TaNiTe2 (part 2)    
ABC3_oP20_53_e_g_hi-001.{ABC}                    oP20      Pmna                -         CuBrSe3             #53       CuBrSe3 (part 2)    
ABC3_oP20_54_e_d_cf-001.{ABC}                    oP20      Pcca                -         BiGaO3              #54       BiGaO3 (part 2)     
ABC_oP6_59_a_a_a-001.{ABC}                       oP6       Pmmn                -         CNCl                #59       CNCl (part 1)       
ABC_oP6_59_a_b_a-001.{ABC}                       oP6       Pmmn                E0_{5}    FeOCl               #59       FeOCl (part 2)
A2B4C_oP28_62_ac_2cd_c-001.{ABC}                 oP28      Pnma                S1_{2}    Mg2SiO4             #62       Forsterite (part 2) 
A2B4C_oP28_62_ac_2cd_c-002.{ABC}                 oP28      Pnma                -         Ag2CrO4             #62       Ag2CrO4 ternary oxide (R. Friedrich), ICSD #16298
AB2C3_oP24_62_c_d_cd-001.{ABC}                   oP24      Pnma                E9_{e}    CuFe2S3             #62       Cubanite (part 2)   
AB2C_oP16_62_c_2c_c-001.{ABC}                    oP16      Pnma                F5_6      CuSbS2              #62       Chalcostibite (part 1)
AB3C_oP20_62_c_cd_a-001.{ABC}                    oP20      Pnma                -         CaTiO3              #62       CaTiO3 Pnma Perovskite (part 1)
AB3C_oP20_62_c_cd_a-002.{ABC}                    oP20      Pnma                -         CaO3Zr              #62       CaO3Zr ternary oxide (R. Friedrich), ICSD #97463
AB4C_oP24_62_c_2cd_c-001.{ABC}                   oP24      Pnma                H0_{2}    BaSO4               #62       Barite (part 2)     
A2BC3_oC24_63_e_c_cg-001.{ABC}                   oC24      Cmcm                -         KFe2S3              #63       Rasvumite (part 2)  
A43B5C17_oC260_63_c8fg6h_cfg_ce3f2h-001.{ABC}    oC260     Cmcm                -         La43Ni17Mg5         #63       La43Ni17Mg5 (part 2)
AB2C_oC16_63_c_2c_c-001.{ABC}                    oC16      Cmcm                -         SrCuO2              #63       SrCuO2 (part 1)     
AB3C_oC20_63_a_cf_c-001.{ABC}                    oC20      Cmcm                -         MgSiO3              #63       Post-perovskite (part 2)
AB4C_oC24_63_a_fg_c-001.{ABC}                    oC24      Cmcm                -         MgO4S               #63       MgSO4 (part 2)      
AB4C_oC24_63_c_fg_c-001.{ABC}                    oC24      Cmcm                H0_{1}    CaSO4               #63       Anhydrite (part 2)  
A2B2C_oC80_64_efg_efg_df-001.{ABC}               oC80      Cmce                -         MgB2C2              #64       MgB2C2 (part 1)     
A2B4C_oC28_66_l_kl_a-001.{ABC}                   oC28      Cccm                -         SrAl2Se4            #66       SrAl2Se4 (part 2)   
A2BC_oC16_67_ag_b_g-001.{ABC}                    oC16      Cmme                -         Al2CuIr             #67       Al2CuIr (part 2)    
ABC2_oC16_67_b_g_ag-001.{ABC}                    oC16      Cmme                -         HoCuP2              #67       HoCuP2 (part 2)     
ABC_oI12_71_h_j_g-001.{ABC}                      oI12      Immm                -         NbPS                #71       NbPS (part 2)       
A2BC_tP16_76_2a_a_a-001.{ABC}                    tP16      P4_{1}              -         LaRhC2              #76       LaRhC2 (part 2)     
A2B7C2_tP88_78_4a_14a_4a-001.{ABC}               tP88      P4_{3}              -         Sr2As2O7            #78       Sr2As2O7 (part 2)   
A2BC2_tI20_79_c_2a_c-001.{ABC}                   tI20      I4                  -         TlZn2Sb2            #79       TlZn2Sb2 (part 2)   
A2BC4_tI14_82_bc_a_g-001.{ABC}                   tI14      I-4                 E3        CdAl2S4             #82       CdAl2S4 (part 1)    
AB4C_tI12_82_c_g_a-001.ABC/BCA/CAB               tI12      I-4                 H0_7      BPO4                #82       BPO4 (part 1)       
ABC_tP24_91_d_d_d-001.{ABC}                      tP24      P4_{1}22            -         ThBC                #91       ThBC (part 2)       
A14B3C5_tP44_94_c3g_ad_bg-001.{ABC}              tP44      P4_{2}2_{1}2        -         Na5Fe3F14           #94       Na5Fe3F14 (part 2)  
A6B2C_tP18_94_eg_c_a-001.{ABC}                   tP18      P4_{2}2_{1}2        -         Li2MoF6             #94       Li2MoF6 (part 2)    
ABC_tP24_95_d_d_d-001.{ABC}                      tP24      P4_{3}22            -         ThBC                #95       ThBC (part 2)       
AB8C2_tI44_97_e_2k_cd-001.{ABC}                  tI44      I422                -         Ta2Se8I             #97       Ta2Se8I (part 2)    
A3BC_tP5_99_bc_a_b-001.{ABC}                     tP5       P4mm                -         Pb(Zr0.52Ti0.48)O3  #99       Tetragonal PZT [PbO3] (part 1)
A3B11C6_tP40_100_ac_bc2d_cd-001.{ABC}            tP40      P4bm                -         Ce3Si6N11           #100      Ce3Si6N11 (part 2)  
A7B7C2_tP32_101_bde_ade_d-001.{ABC}              tP32      P4_{2}cm            -         gamma-MgNiSn        #101      gamma-MgNiSn (part 2)
A5B5C4_tP28_104_ac_ac_c-001.{ABC}                tP28      P4nc                -         Ba5In4Bi5           #104      Ba5In4Bi5 (part 2)  
AB6C4_tP22_104_a_2ac_c-001.{ABC}                 tP22      P4nc                -         Tl4HgI6             #104      Tl4HgI6 (part 2)    
A2BC2_tP20_105_f_ac_2e-001.{ABC}                 tP20      P4_{2}mc            -         BaGe2As2            #105      BaGe2As2 (part 2)   
ABC_tI12_109_a_a_a-001.ABC/BCA/CAB               tI12      I4_{1}md            -         LaPtSi              #109      LaPtSi (part 2)     
A2BC8_tI176_110_2b_b_8b-001.{ABC}                tI176     I4_{1}cd            -         Be[BH4]2            #110      Be2 (part 2)        
AB4C_tP12_112_b_n_e-001.{ABC}                    tP12      P-42c               -         alpha-CuAlCl4       #112      alpha-CuAlCl4 (part 2)
A4BC2_tI28_120_i_d_e-001.{ABC}                   tI28      I-4c2               -         KAu4Sn2             #120      KAu4Sn2 (part 2)    
ABC2_tI16_122_a_b_d-001.{ABC}                    tI16      I-42d               E1_1      CuFeS2              #122      Chalcopyrite (part 1)
AB5C_tP7_123_b_ci_a-001.{ABC}                    tP7       P4/mmm              -         HoCoGa5             #123      HoCoGa5 (part 1)    
ABC2_tP4_123_d_a_f-001.{ABC}                     tP4       P4/mmm              -         CaCuO2              #123      CaCuO2 (part 1)     
AB4C_tP12_124_a_m_c-001.ABC/BCA/CAB              tP12      P4/mcc              -         Nb4CoSi             #124      Nb4CoSi (part 2)    
ABC4_tP12_125_a_b_m-001.{ABC}                    tP12      P4/nbm              -         KCeSe4              #125      KCeSe4 (part 2)     
A2BC4_tP28_126_cd_e_k-001.{ABC}                  tP28      P4/nnc              -         BiAl2S4             #126      BiAl2S4 (part 2)    
A6B2C_tP18_128_eh_d_a-001.{ABC}                  tP18      P4/mnc              -         K2SnCl6             #128      K2SnCl6 (part 2)    
A7B2C_tP40_128_egi_h_e-001.{ABC}                 tP40      P4/mnc              E9_{a}    FeCu2Al7            #128      FeCu2Al7 (part 2)   
ABC_tP6_129_c_a_c-001.{ABC}                      tP6       P4/nmm              E0_1      PbFCl               #129      Matlockite (part 1) 
A2BC4_tP28_130_f_c_g-001.{ABC}                   tP28      P4/ncc              -         CuBi2O4             #130      CuBi2O4 (part 2)    
AB6C_tP16_132_d_io_a-001.{ABC}                   tP16      P4_{2}/mcm          -         AgUF6               #132      AgUF6 (part 2)      
A4B2C_tP28_135_gh_h_d-001.{ABC}                  tP28      P4_{2}/mbc          -         ZnSb2O4             #135      ZnSb2O4 (part 2)    
A4BC4_tP18_137_g_b_g-001.{ABC}                   tP18      P4_{2}/nmc          -         CeCo4B4             #137      CeCo4B4 (part 2)    
AB2C4_tI14_139_a_e_ce-001.{ABC}                  tI14      I4/mmm              -         (La,Ba)2CuO4        #139      0201 [(La,Ba)2CuO4] (part 1)
AB2C4_tI14_139_a_e_ce-002.{ABC}                  tI14      I4/mmm              -         O4Sr2Ti             #139      O4Sr2Ti ternary oxide (R. Friedrich), ICSD #157402 
A4BC_tI24_141_h_b_a-001.{ABC}                    tI24      I4_{1}/amd          -         ZrSiO4              #141      Zircon (part 1)     
ABC4_tI96_142_e_ab_2g-001.{ABC}                  tI96      I4_{1}/acd          -         PPrS4               #142      PPrS4 (part 1)      
A4B14C3_hP21_143_bd_ac4d_d-001.{ABC}             hP21      P3                  -         Ta3Al4O13[OH]       #143      Simpsonite (part 2) 
A4B6C_hP11_143_bd_2d_a-001.{ABC}                 hP11      P3                  -         ScRh6P4             #143      ScRh6P4 (part 2)    
A3BC_hR5_146_b_a_a-001.{ABC}                     hR5       R3                  -         gamma-Ag3SI         #146      gamma-Ag3SI (part 2)
ABC3_hR10_146_2a_2a_2b-001.{ABC}                 hR10      R3                  -         FePSe3              #146      FePSe3 (part 2)     
A2B4C_hR42_148_2f_4f_f-001.{ABC}                 hR42      R-3                 S1_{3}    Be2SiO4             #148      Phenakite (part 2)  
AB3C_hR10_148_c_f_c-001.ABC/BCA/CAB              hR10      R-3                 -         FeTiO3              #148      Ilmenite (part 1)   
A5B6C2_hP13_157_2ac_2c_b-001.{ABC}               hP13      P31m                -         Ag5Pb2O6            #157      Ag5Pb2O6 (part 2)   
ABC_hR3_160_a_a_a-001.ABC/BAC/CAB                hR3       R3m                 F0_{2}    COS                 #160      Carbonyl Sulphide (part 2)
ABC3_hR10_161_a_a_b-001.{ABC}                    hR10      R3c                 -         LiNbO3              #161      Ferroelectric LiNbO3 (part 1)
ABC2_hP4_164_a_b_d-001.{ABC/BCA/CAB              hP4       P-3m1               -         CuNiSb2             #164      CuNiSb2 (part 2)    
ABC2_hR4_166_a_b_c-001.{ABC}                     hR4       R-3m                F5_1      CrNiS2              #166      Caswellsilverite (part 1)
ABC2_hR4_166_a_b_c-002.ABC/BCA/CAB               hR4       R-3m                -         ABC2                #166      rocksalt derivative (Y. Lederer)
ABC2_hR4_166_a_b_c-003.ABC/BCA/CAB               hR4       R-3m                -         AlLiO2              #166      AlLiO2 ternary oxide (R. Friedrich), ICSD #28288 
ABC2_hR4_166_a_b_c-004.{ABC}                     hR4       R-3m                -         CuFeO2              #166      CuFeO2 ternary oxide (R. Friedrich), ICSD #246912 
ABC_hR6_166_c_c_c-001.{ABC}                      hR6       R-3m                -         SmSI                #166      SmSI (part 2)       
AB3C_hR10_167_b_e_a-001.{ABC}                    hR10      R-3c                -         PrNiO3              #167      PrNiO3 (part 2)     
ABC2_hR24_167_e_e_2e-001.{ABC}                   hR24      R-3c                F5_{13}   KBO2                #167      KBO2 (part 2)       
ABC3_hR10_167_a_b_e-001.{ABC}                    hR10      R-3c                -         LiNbO3              #167      Paraelectric LiNbO3 (part 1)
ABC3_hR10_167_a_b_e-002.{ABC}                    hR10      R-3c                G0_1      CaCO3               #167      Calcite (part 1)    
A2B13C4_hP57_168_d_c6d_2d-001.{ABC}              hP57      P6                  -         K2Ta4O9F4           #168      K2Ta4O9F4 (part 2)  
AB4C_hP72_168_2d_8d_2d-001.{ABC}                 hP72      P6                  -         Al[PO4]             #168      Al (part 2)         
A10B2C_hP39_171_5c_c_a-001.{ABC}                 hP39      P6_{2}              -         Sr[S2O6][H2O]4      #171      Sr4 (part 2)        
A10B2C_hP39_172_5c_c_a-001.{ABC}                 hP39      P6_{4}              -         Sr[S2O6][H2O]4      #172      Sr4 (part 2)        
A12B7C2_hP21_174_2j2k_ajk_cf-001.{ABC}           hP21      P-6                 -         Fe12Zr2P7           #174      Fe12Zr2P7 (part 2)  
ABC_hP12_174_cj_fk_aj-001.ABC/BCA/CAB            hP12      P-6                 -         GdSI                #174      GdSI (part 2)       
A8B7C6_hP21_175_ck_aj_k-001.{ABC}                hP21      P6/m                -         Nb7Ru6B8            #175      Nb7Ru6B8 (part 2)   
ABC_hP36_175_jk_jk_jk-001.{ABC}                  hP36      P6/m                -         Mg[NH]              #175      Mg (part 2)         
A3B3C_hP14_176_h_h_c-001.{ABC}                   hP14      P6_{3}/m            -         Fe3Te3Tl            #176      Fe3Te3Tl (part 3)
ABC_hP3_183_a_a_a-001.ABC/BCA/CAB                hP3       P6mm                -         AuCN                #183      AuCN (part 2)       
AB4C_hP72_184_d_4d_d-001.ABC/BCA/CAB             hP72      P6cc                -         Al[PO4]             #184      Al (part 2)         
A3BC_hP30_185_cd_c_ab-001.{ABC}                  hP30      P6_{3}cm            -         KNiCl3              #185      KNiCl3 (part 2)     
A5B3C_hP18_186_2a3b_2ab_b-001.{ABC}              hP18      P6_{3}mc            E9_4      Al5C3N              #186      Al5C3N (part 1)     
ABC_hP3_187_a_d_f-001.ABC/BCA/CAB                hP3       P-6m2               -         BaPtSb              #187      BaPtSb (part 1)     
A3BC_hP10_188_k_c_a-001.ABC/BCA/CAB              hP10      P-6c2               -         LiScI3              #188      LiScI3 (part 2)
AB9C4_hP28_188_e_kl_ak-001.{ABC}                 hP28      P-6c2               S3_{2}    BaSi4O9             #188      BaSi4O9 (part 2)
AB4C_hP6_191_a_h_b-001.ABC/BCA/CAB               hP6       P6/mmm              -         AlB4Mg              #191      AlB4Mg (part 1)     
A9B3C_hP26_194_hk_h_a-001.{ABC}                  hP26      P6_{3}/mmc          E9_{c}    Al9Mn3Si            #194      Al9Mn3Si (part 2)   
AB3C4_hP16_194_c_af_ef-001.{ABC}                 hP16      P6_{3}/mmc          -         AlN3Ti4             #194      AlN3Ti4 (part 1)    
ABC2_hP8_194_d_a_f-001.{ABC}                     hP8       P6_{3}/mmc          -         AlCCr2              #194      AlCCr2 (part 1)     
ABC_hP6_194_c_d_a-001.ABC/BCA/CAB                hP6       P6_{3}/mmc          -         LiBC                #194      LiBC (part 1)       
A12BC4_cP34_195_2j_ab_2e-001.{ABC}               cP34      P23                 -         PrRu4P12            #195      PrRu4P12 (part 2)   
A12B2C_cF60_196_h_bc_a-001.{ABC}                 cF60      F23                 -         Cu2Fe[CN]6          #196      Cu2Fe6 (part 2)     
ABC3_cP20_198_a_a_b-001.{ABC}                    cP20      P2_{1}3             G3        NaClO3              #198      Sodium Chlorate (part 2)
ABC_cP12_198_a_a_a-001.{ABC}                     cP12      P2_{1}3             F0_1      NiSSb               #198      Ullmanite (part 1)  
AB3C_cP60_201_be_fh_g-001.{ABC}                  cP60      Pn-3                -         KSbO3               #201      KSbO3 (part 2)
A6B6C_cF104_202_h_h_c-001.{ABC}                  cF104     Fm-3                -         KB6H6               #202      KB6H6 (part 2)      
AB3C16_cF160_203_a_bc_eg-001.{ABC}               cF160     Fd-3                -         Rb3AsSe16           #203      Rb3AsSe16 (part 2)
AB32C48_cI162_204_a_2efg_2gh-001.{ABC}           cI162     Im-3                -         Mg32(Al,Zn)49       #204      Bergman [Mg32(Al,Zn)49] (part 1)
A2B3C6_cP264_205_2d_ab2c2d_6d-001.{ABC}          cP264     Pa-3                -         Ca3Al2O6            #205      Ca3Al2O6 (part 2)   
AB3C2_cI96_206_c_e_ad-001.{ABC}                  cI96      Ia-3                E9_{d}    AlLi3N2             #206      AlLi3N2 (part 2)    
AB3C6_cI80_206_a_d_e-001.{ABC}                   cI80      Ia-3                D5_3      (Mn,Fe)2O3          #206      Bixbyite (part 1)   
A24BC_cF104_209_j_a_b-001.{ABC}                  cF104     F432                -         F6KP                #209      F6KP (part 2)       
A12B6C_cF608_210_4h_2h_e-001.{ABC}               cF608     F4_{1}32            -         Te[OH]6             #210      Te6 (part 2)        
A3B3C_cI56_214_g_h_a-001.{ABC}                   cI56      I4_{1}32            -         Ca3PI3              #214      Ca3PI3 (part 2)     
A3BC2_cI48_214_f_a_e-001.{ABC}                   cI48      I4_{1}32            -         Ag3AuTe2            #214      Petzite (part 2)    
A3B4C_cP8_215_d_e_a-001.{ABC}                    cP8       P-43m               H2_4      Cu3S4V              #215      Sulvanite (part 1)  
AB3C4_cP8_215_a_c_e-001.{ABC}                    cP8       P-43m               -         Cu3AsS4             #215      Cubic Lazarevicite (part 1)
ABC_cF12_216_b_c_a.ABC/BCA/CAB                   cF12      F-43m               C1_b      AgAsMg              #216      Half-Heusler (part 1)
A3B4C_cP16_218_c_e_a-001.{ABC}                   cP16      P-43n               H2_{1}    Ag3[PO4]            #218      Ag3 (part 2)
A2B3C6_cP33_221_cd_ag_fh-001.{ABC}               cP33      Pm-3m               E9_{1}    Ca3Al2O6            #221      Ca3Al2O6 (part 2)   
AB3C_cP5_221_a_c_b.{ABC}                         cP5       Pm-3m               E2_1      CaTiO3              #221      Cubic Perovskite (part 1)
A5B3C16_cP96_222_ce_d_fi-001.{ABC}               cP96      Pn-3n               -         Ce5Mo3O16           #222      Ce5Mo3O16 (part 2)  
A6B2C_cF36_225_e_c_a-001.{ABC}                   cF36      Fm-3m               J1_{1}    K2PtCl6             #225      K2PtCl6 (part 2)    
A9B16C7_cF128_225_acd_2f_be-001.{ABC}            cF128     Fm-3m               -         Cr9Fe16Ni7          #225      Model of Ferrite (part 1)
AB18C8_cF108_225_a_eh_f-001.{ABC}                cF108     Fm-3m               -         CrFe18Ni8           #225      Model of Austenite (part 1)
AB2C_cF16_225_a_c_b.ABC/BCA/CAB                  cF16      Fm-3m               L2_1      AlCu2Mn             #225      Heusler (part 1), B2 derivative (Y. Lederer), equivalent to T0001 [HTQC] 
A2B2C7_cF88_227_c_d_af-001.{ABC}                 cF88      Fd-3m               -         Eu2Ir2O7            #227      Pyrochlore Iridate (part 2)
A2BC4_cF56_227_d_a_e-001.{ABC}                   cF56      Fd-3m               H1_1      Al2MgO4             #227      Spinel (part 1)     
AB3C3_cF112_227_c_de_f-001.{ABC}                 cF112     Fd-3m               E9_3      Fe3W3C              #227      Fe3W3C (part 1)     
AB12C3_cI32_229_a_h_b-001.{ABC}                  cI32      Im-3m               -         CrFe12Ni3           #229      Model of Austenite (part 1)
AB4C3_cI16_229_a_c_b.{ABC}                       cI16      Im-3m               -         CrFe4Ni3            #229      Model of Ferrite (part 1)

*****************************************************************************************************************************************************
ANRL - QUATERNARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes (part)
----------------------------------------------------------------------------------------------------------------------------------------------------
ABC6D2_mC40_15_e_e_3f_f-001.{ABCD}               mC40      C2/c                -         CaFeO6Si2           #15       Esseneite (part 1)  
AB7CD2_oI44_24_a_b3d_c_ac-001.{ABCD}             oI44      I2_{1}2_{1}2_{1}    -         Na2MgAlF7           #24       Weberite (part 2)   
A6B4C16D_oP108_27_abcd4e_4e_16e_e-001.{ABCD}     oP108     Pcc2                -         Ca4Al6O16S          #27       Ca4Al6O16S (part 2) 
A2B3C7D_oP13_47_t_aq_eqrs_h-001.{ABCD}           oP13      Pmmm                -         YBa2Cu3O7-x         #47       1212C [YBa2Cu3O7-x] (part 1)
A2B8CD_oP24_48_k_2m_d_b-001.{ABCD}               oP24      Pnnn                -         alpha-RbPr[MoO4]2   #48       alpha-RbPr2 (part 2)
AB2C8D_oP24_49_g_q_2qr_e-001.{ABCD}              oP24      Pccm                -         CsPr[MoO4]2         #49       CsPr2 (part 2)      
ABCD_oP16_57_d_c_d_d-001.{ABCD}                  oP16      Pbcm                F5_9      KCNS                #57       KCNS (part 1)       
ABCD3_oI48_73_d_e_e_ef-001.{ABCD}                oI48      Ibca                -         KAg[CO3]            #73       KAg (part 2)        
AB2C12D4_tP76_75_2a2b_2d_12d_4d-001.{ABCD}       tP76      P4                  -         BaCr2Ru4O12         #75       BaCr2Ru4O12 (part 2)
A2B6CD7_tP64_77_2d_6d_d_ab6d-001.{ABCD}          tP64      P4_{2}              -         MgB2O(OH)6          #77       Pinnoite (part 2)   
A17BC4D_tP184_89_17p_p_4p_io-001.{ABCD}          tP184     P422                -         C17FeO4Pt           #89       C17FeO4Pt (part 2)  
A4B2C13D_tP40_90_g_d_cef2g_c-001.{ABCD}          tP40      P42_{1}2            -         Na4Ti2Si8O22[H2O]4  #90       Na4Ti2Si8O224 (part 2)
A2B8CD_tI24_97_d_k_a_b-001.ABCD/ACDB/ADBC/       tI24      I422                -         NaGdCu2F8           #97       NaGdCu2F8 (part 2)  
BADC/BCAD/BDCA/CABD/CBDA/CDAB/DACB/DBAC/DCBA
A2B8C2D_tP26_100_c_abcd_c_a-001.{ABCD}           tP26      P4bm                -         Ba2TiSi2O8          #100      Fresnoite (part 2)  
A3BC3D_tP64_106_3c_c_3c_c-001.{ABCD}             tP64      P4_{2}bc            -         NaZn[OH]3           #106      NaZn3 (part 2)      
A2BC7D2_tP24_113_e_a_cef_e-001.{ABCD}            tP24      P-42_{1}m           S5_{3}    Ca2MgSi2O7          #113      Akermanite (part 2) 
A2BC4D_tI16_121_d_a_i_b-001.ABCD/ACDB/ADBC/      tI16      I-42m               H2_6      Cu2FeS4Sn           #121      Stannite (part 1)   
BADC/BCAD/BDCA/CABD/CBDA/CDAB/DACB/DBAC/DCBA
A4BC4D_tP10_123_gh_a_i_d-001.{ABCD}              tP10      P4/mmm              -         CaRbFe4As4          #123      CaRbFe4As4 (part 2)
ABCD_tP8_129_c_b_a_c-001.{ABCD}                  tP8       P4/nmm              -         AsCuSiZr            #129      AsCuSiZr (part 1)   
A2B2C4D_tP18_132_e_i_o_d-001.{ABCD}              tP18      P4_{2}/mcm          -         Rb2TiCu2Se4         #132      Rb2TiCu2S4 (part 2) 
AB4C7D_hP26_159_b_ac_a2c_b-001.{ABCD}            hP26      P31c                -         YbBaCo4O7           #159      YbBaCo4O7 (part 2)  
AB2CD2_hP36_163_h_i_bf_i-001.{ABCD}              hP36      P-31c               F5_10     KAg(CN)2            #163      KAg2 (part 1)       
A8BC3D6_hP18_189_bfh_a_g_i-001.{ABCD}            hP18      P-62m               E9_{b}    pi-FeMg3Al8Si6      #189      pi-FeMg3Al8Si6 (part 2)
A9BC3D5_hP18_189_fi_a_g_bh-001.{ABCD}            hP18      P-62m               -         pi-FeMg3Al9Si5      #189      pi-FeMg3Al9Si5 (part 2)
A2B3C18D6_hP58_192_c_f_lm_l-001.{ABCD}           hP58      P6/mcc              G3_{1}    Be3Al2Si6O18        #192      Beryl (part 2)      
A12B36CD12_cF488_210_h_3h_a_fg-001.{ABCD}        cF488     F4_{1}32            -         MgB12H12[H2O]12     #210      MgB12H1212 (part 2)
ABCD_cF16_216_c_d_b_a.ABCD/BCAD/CABD             cF16      F-43m               -         LiMgAuSn            #216      Quaternary Heusler (part 2)
A7BC3D13_cF192_219_de_b_c_ah-001.{ABCD}          cF192     F-43c               -         Mg3B7ClO13          #219      Boracite (part 2)   
AB11CD3_cP16_221_a_dg_b_c-001.{ABCD}             cP16      Pm-3m               -         CrFe11MoNi3         #221      Model of Ferrite (part 1)
AB27CD3_cP32_221_a_dij_b_c-001.{ABCD}            cP32      Pm-3m               -         CrFe2525Ni6         #221      Model of Austenite (part 1)
A5BCD6_cF416_228_eg_c_b_h-001.{ABCD}             cF416     Fd-3c               -         CuCrCl5[NH3]6       #228      CuCrCl56 (part 2)   
A2B3C12D3_cI160_230_a_c_h_d-001.{ABCD}           cI160     Ia-3d               S1_{4}    Co3Al2Si3O12        #230      Garnet (part 2)     

*****************************************************************************************************************************************************
ANRL - QUINARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes (part)
----------------------------------------------------------------------------------------------------------------------------------------------------
A8B2C12D2E_oI50_23_bcfk_i_3k_j_a-001.{ABCDE}     oI50      I222                -         Cu8(Fe,Zn)3Sn2S12   #23       Stannoidite (part 2)
AB4C17D4E_tP54_90_a_g_c4g_g_c-001.{ABCDE}        tP54      P42_{1}2            -         BaCu4[VO][PO4]4     #90       BaCu44 (part 2)     
AB32CD4E8_tP184_93_i_16p_af_2p_4p-001.{ABCDE}    tP184     P4_{2}22            -         AsPh4CeS8P4Me8      #93       AsPh4CeS8P4Me8 (part 2)
A2B3C3DE7_hP48_145_2a_3a_3a_a_7a-001.{ABCDE}     hP48      P3_{2}              -         NaCa3[CO3]2F3[H2O]  #145      Sheldrickite (part 2)
A2BCD3E6_cF208_203_e_c_d_f_g-001.{ABCDE}         cF208     Fd-3                -         Na3Co(CO3)2Cl       #203      Pyrochlore (part 2) 
A4B2C6D16E_cF232_203_e_d_f_eg_a-001.{ABCDE}      cF232     Fd-3                -         Na6Mg2(SO4)(CO3)4   #203      Tychite (part 2)    
A6B2CD6E_cP64_208_m_ad_b_m_c-001.{ABCDE}         cP64      P4_{2}32            -         Cs2ZnFe[CN]6        #208      Cs2ZnFe6 (part 2)   

*****************************************************************************************************************************************************
ANRL - BCC-SQS

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
AB_aP16_2_4i_4i-001.AB                           aP16      P-1                 -         AB                  #2        bcc-SQS [10.1103/PhysRevB.69.214202] (O. Levy)
A5B11_mP16_6_2abc_2a3b3c-001.{AB}                mP16      Pm                  -         A5B11               #6        bcc-SQS [10.1103/PhysRevB.94.224104] (O. Levy)      
AB3_mC32_8_4a_12a-001.{AB}                       mC32      Cm                  -         AB3                 #8        bcc-SQS [10.1103/PhysRevB.69.214202] (O. Levy)
AB3_mC32_8_4a_4a4b-001.{AB}                      mC32      Cm                  -         AB3                 #8        bcc-SQS [10.1103/PhysRevB.94.224104] (O. Levy)
A3B13_oC32_38_ac_a2bcdef-001.{AB}                oC32      Amm2                -         A3B13               #38       bcc-SQS [10.1103/PhysRevB.94.224104] (O. Levy)
A3B5_oC32_38_abce_abcdf-001.{AB}                 oC32      Amm2                -         A3B5                #38       bcc-SQS [10.1103/PhysRevB.94.224104] (O. Levy)
AB7_hR16_166_c_c2h-001.{AB}                      hR16      R-3m                -         AB7                 #166      bcc-SQS [10.1103/PhysRevB.94.224104] (O. Levy)

*****************************************************************************************************************************************************
ANRL - miscellaneous

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
A2BCD4_tI16_82_ac_b_d_g-001.ABCD/BCAD/           tI16      I-4                 -         Cu2(Zn,Fe)SnS4      #82       Kesterite, 12 protos [Can. Mineral., 16, 131-137 (1978)] (C. Oses)
CABD/DBAC/BADC/ADBC/ACDB/CDAB/DACB/
DCBA/CBDA/BDCA   
A2B3_cI80_206_ad_e.{AB}                          tI80      Ia-3                -         NZn                 #206      NZn (ICSD #84918) (R. Friedrich)

*****************************************************************************************************************************************************
ANRL - fcc/sc derivative structures - BINARIES (Y. Lederer)

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
AB3_mC8_12_a_di-001.{AB}                         mC8       C2/m                -         AB3                 #12       fcc derivative, equivalent to 16,18 [HTQC] 
AB_mC8_12_i_i-001.AB                             mC8       C2/m                -         AB                  #12       fcc derivative, equivalent to 17 [HTQC]
AB3_oP4_47_a_ct-001.{AB}                         oP4       Pmmm                -         AB3                 #47       fcc derivative, equivalent to 19,21 [HTQC]
A3B_oP4_47_cr_a-001.{AB}                         oP4       Pmmm                -         A3B                 #47       sc derivative
AB_oP4_51_e_e-001.AB                             oP4       Pmma                -         AB                  #51       sc derivative
AB_oP4_59_a_a-001.AB                             oP4       Pmmn                -         AB                  #59       fcc derivative, equivalent to 20 [HTQC]
AB2_oC6_65_a_i-001.{AB}                          oC6       Cmmm                -         AB2                 #65       sc derivative
A3B_oC8_65_ai_b-001.{AB}                         oC8       Cmmm                -         A3B                 #65       sc derivative
A3B_oC8_65_bj_a-001.{AB}                         oC8       Cmmm                -         A3B                 #65       sc derivative
AB_oC8_65_i_i-001.AB                             oC8       Cmmm                -         AB                  #65       sc derivative
AB_tP2_123_a_b-001.AB                            tP2       P4/mmm              -         AB                  #123      sc derivative
AB_tP2_123_a_c-001.AB                            tP2       P4/mmm              -         AB                  #123      sc derivative 
A2B_tP3_123_g_a-001.{AB}                         tP3       P4/mmm              -         A2B                 #123      sc derivative 
A3B_tP4_123_abc_d-001.{AB}                       tP4       P4/mmm              -         A3B                 #123      sc derivative 
A3B_tP4_123_ag_b-001.{AB}                        tP4       P4/mmm              -         A3B                 #123      sc derivative
A3B_tP4_123_cf_a-001.{AB}                        tP4       P4/mmm              -         A3B                 #123      sc derivative 
AB3_tP4_123_a_bh-001.{AB}                        tP4       P4/mmm              -         AB3                 #123      fcc derivative, equivalent to 13,15 [HTQC]
AB_tP4_123_g_g-001.AB                            tP4       P4/mmm              -         AB                  #123      sc derivative
A3B_tI8_139_ae_b-001.{AB}                        tI8       I4/mmm              -         A3B                 #139      sc derivative 
A2B_hP3_164_d_a-001.{AB}                         hP3       P-3m1               -         A2B                 #164      fcc derivative, equivalent to 9,10 [HTQC]
A2B_hP3_164_d_a-002.{AB}                         hP3       P-3m1               -         A2B                 #164      sc derivative 
A3B_hR4_166_bc_a-001.{AB}                        hR4       R-3m                -         A3B                 #166      sc derivative 
AB3_hR4_166_a_bc-001.{AB}                        hR4       R-3m                -         AB3                 #166      fcc derivative, equivalent to 27,29 [HTQC]
AB_hR4_166_c_c-001.AB                            hR4       R-3m                -         AB                  #166      fcc derivative, equivalent to 28 [HTQC]
AB_hR4_166_c_c-002.AB                            hR4       R-3m                -         AB                  #166      sc derivative

*****************************************************************************************************************************************************
ANRL - rocksalt/B2 derivative structures - TERNARIES (Y. Lederer)

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
AB3C4_mC16_12_a_di_2i-001.{ABC}                  mC16      C2/m                -         AB3C4               #12       rocksalt derivative, AB3_mC8_12_a_di-001 with C "in between"
ABC2_mC16_12_i_i_adi-001.ABC/BCA/CAB             mC16      C2/m                -         ABC2                #12       rocksalt derivative, AB_mC8_12_i_i-001 with C "in between"
AB3C4_oP8_47_a_ct_egs-001.{ABC}                  oP8       Pmmm                -         AB3C4               #47       rocksalt derivative, AB3_oP4_47_a_ct-001 with C "in between"
A3BC4_oP8_47_eq_g_bdt-001.{ABC}                  oP8       Pmmm                -         A3BC4               #47       B2 derivative, A3B_oP4_47_cr_a-001 with C "in between", shifted origin
ABC2_oP8_51_e_e_2f-001.ABC/BCA/CAB               oP8       Pmma                -         ABC2                #51       B2 derivative, AB_oP4_51_e_e-001 with C "in between"
ABC2_oP8_59_a_a_2b-001.ABC/BCA/CAB               oP8       Pmmn                -         ABC2                #59       rocksalt derivative, AB_oP4_59_a_a-001 with C "in between"
ABC2_oC16_63_c_c_g-001.ABC/BCA/CAB               oC16      Cmcm                -         ABC2                #63       B2 derivative, AB_oC8_63_c_c-002 with C "in between" 
AB2C3_oC12_65_a_i_cj-001.{ABC}                   oC12      Cmmm                -         AB2C3               #65       B2 derivative, AB2_oC6_65_a_i-001 with C "in between" 
A3BC4_oC16_65_ai_b_q-001.{ABC}                   oC16      Cmmm                -         A3BC4               #65       B2 derivative, A3B_oC8_65_ai_b-001 with C "in between"
A3BC4_oC16_65_bj_a_eh-001.{ABC}                  oC16      Cmmm                -         A3BC4               #65       B2 derivative, A3B_oC8_65_bj_a-001 with C "in between"
AB3C4_oC16_65_a_bf_hi-001.{ABC}                  oC16      Cmmm                -         AB3C4               #65       rocksalt derivative, AB3_oC8_65_a_bf-001 with C "in between" 
ABC2_oC16_65_i_i_fh-001.ABC/BCA/CAB              oC16      Cmmm                -         ABC2                #65       B2 derivative, AB_oC8_65_i_i-001 with C "in between"
AB2C3_oI12_71_a_e_df-001.{ABC}                   oI12      Immm                -         AB2C3               #71       rocksalt derivative, AB2_oI6_71_a_g-001 with C "in between", shifted origin 
ABC2_tP4_123_a_b_h-001.ABC/BCA/CAB               tP4       P4/mmm              -         ABC2                #123      B2 derivative, AB_tP2_123_a_b with C "in between"; equivalent to TBCC011 [HTQC]
ABC2_tP4_123_a_c_e-001.ABC/BCA/CAB               tP4       P4/mmm              -         ABC2                #123      B2 derivative, AB_tP2_123_a_c-001 with C "in between"; equivalent to TBCC016 [HTQC]
ABC2_tP4_123_a_d_bc-001.ABC/BCA/CAB              tP4       P4/mmm              -         ABC2                #123      rocksalt derivative, AB_tP2_123_a_d-001 with C "in between" 
A2BC3_tP6_123_g_b_ch-001.{ABC}                   tP6       P4/mmm              -         A2BC3               #123      B2 derivative, A2B_tP3_123_g_a-001 with C "in between", shifted origin
A3BC4_tP8_123_abc_d_i-001.{ABC}                  tP8       P4/mmm              -         A3BC4               #123      B2 derivative, A3B_tP4_123_abc_d-001 with C "in between"
A3BC4_tP8_123_ag_b_2h-001.{ABC}                  tP8       P4/mmm              -         A3BC4               #123      B2 derivative, A3B_tP4_123_ag_b-001 with C "in between" 
A3BC4_tP8_123_cf_a_k-001.{ABC}                   tP8       P4/mmm              -         A3BC4               #123      B2 derivative, A3B_tP4_123_cf_a-001 with C "in between"
AB3C4_tP8_123_a_bh_cdg-001.{ABC}                 tP8       P4/mmm              -         AB3C4               #123      rocksalt derivative, AB3_tP4_123_a_bh-001 with C "in between"
ABC2_tP8_123_h_h_abg-001.ABC/BCA/CAB             tP8       P4/mmm              -         ABC2                #123      B2 derivative, AB_tP4_123_g_g-001 with C "in between", shifted origin
ABC2_tP8_129_c_c_2c-001.ABC/BCA/CAB              tP8       P4/nmm              -         ABC2                #129      rocksalt derivative, AB_tP4_129_c_c-002 with C "in between" 
AB2C3_tI12_139_a_e_be-001.{ABC}                  tI12      I4/mmm              -         AB2C3               #139      rocksalt derivative, AB_tI8_139_e_e-002 with C "in between"
A3BC4_tI16_139_ae_b_g-001.{ABC}                  tI16      I4/mmm              -         A3BC4               #139      B2 derivative, A3B_tI8_139_ae_b-001 with C "in between" 
AB3C4_tI16_139_a_bd_ce-001.{ABC}                 tI16      I4/mmm              -         AB3C4               #139      rocksalt derivative, A3B_tI8_139_bd_a-002.BA (swap A and B) with C "in between" 
ABC2_tI16_139_e_e_cd-001.ABC/BCA/CAB             tI16      I4/mmm              -         ABC2                #139      B2 derivative, AB_tI8_139_e_e-002 with C "in between" 
ABC2_tI16_141_a_b_e-001.ABC/BCA/CAB              tI16      I4_{1}/amd          -         ABC2                #141      rocksalt derivative, AB_tI8_141_a_b-002 with C "in between"
A2BC3_hP6_164_d_a_bd-001.{ABC}                   hP6       P-3m1               -         A2BC3               #164      rocksalt derivative, A2B_hP3_164_d_a-001 with C "in between"
A2BC3_hP6_164_d_b_ad-001.{ABC}                   hP6       P-3m1               -         A2BC3               #164      B2 derivative, A2B_hP3_164_d_a-002 with C "in between", shifted origin 
A3BC4_hR8_166_bc_a_2c-001.{ABC}                  hR8       R-3m                -         A3BC4               #166      B2 derivative, A3B_hR4_166_bc_a-001 with C "in between"
AB3C4_hR8_166_a_bc_2c-001.{ABC}                  hR8       R-3m                -         AB3C4               #166      rocksalt derivative, AB3_hR4_166_a_bc-001 with C "in between"
ABC2_hR8_166_c_c_abc-001.ABC/BCA/CAB             hR8       R-3m                -         ABC2                #166      rocksalt derivative, AB_hR4_166_c_c-001 with C "in between"
ABC2_hR8_166_c_c_abc-002.ABC/BCA/CAB             hR8       R-3m                -         ABC2                #166      B2 derivative, AB_hR4_166_c_c-002 with C "in between"
AB3C4_cP8_221_a_c_bd.{ABC}                       cP8       Pm-3m               -         AB3C4               #221      rocksalt derivative, AB3_cP4_221_a_c with C "in between" 

*****************************************************************************************************************************************************
ANRL - oxides - BINARIES (R. Friedrich)

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
AB_mC20_12_b2i_c2i-001.{AB}                      mC20      C2/m                -         OTi                 #12       OTi, ICSD #56694
A2B3_mC20_12_2i_3i-001.{AB}                      mC20      C2/m                -         Ga2O3               #12       Ga2O3, ICSD #83645
A5B3_mC32_12_5i_3i-001.{AB}                      mC32      C2/m                -         O5Ti3               #12       O5Ti3, ICSD #75193
A3B_mP16_14_3e_e-001.{AB}                        mP16      P2_{1}/c            -         O3W                 #14       O3W, ICSD #84848
A2B3_mP20_14_2e_3e-001.{AB}                      mP20      P2_{1}/c            -         Bi2O3               #14       Bi2O3, ICSD #168806
AB_mC8_15_a_e-001.{AB}                           mC8       C2/c                -         CuO                 #15       CuO, ICSD #92368
A4B_mC20_15_2f_e-001.{AB}                        mC20      C2/c                -         O4Os                #15       O4Os, ICSD #23803
A2B5_oP28_19_2a_5a-001.{AB}                      oP28      P2_{1}2_{1}2_{1}    -         As2O5               #19       As2O5, ICSD #654040
A2B_oP24_33_4a_2a-001.{AB}                       oP24      Pna2_{1}            -         O2Sb                #33       O2Sb, ICSD #4109
AB3_oC16_40_b_a2b-001.{AB}                       oC16      Ama2                -         CrO3                #40       CrO3, ICSD #16031
A5B2_oP14_59_a2e_e-001.{AB}                      oP14      Pmmn                -         O5V2                #59       O5V2, ICSD #60767
AB_oC16_64_e_f-001.{AB}                          oC16      Cmce                -         KO                  #64       KO, ICSD #180559
A4B3_tP28_135_gh_dh-001.{AB}                     tP28      P4_{2}/mbc          -         O4Pb3               #135      O4Pb3, ICSD #29094
A2B3_hP15_144_2a_3a-001.{AB}                     hP15      P3_{1}              -         B2O3                #144      B2O3, ICSD #24649
A2B_hR3_166_c_a-001.{AB}                         hR3       R-3m                -         Cs2O                #166      Cs2O, ICSD #27919
AB2_hR6_166_c_2c-001.{AB}                        hR6       R-3m                -         OTl2                #166      OTl2, ICSD #16220
AB_hP12_189_fg_eh-001.{AB}                       hP12      P-62m               -         NaO                 #189      NaO, ICSD #25526
AB_hP8_194_ac_f-001.{AB}                         hP8       P6_{3}/mmc          -         LiO                 #194      LiO, ICSD #152183
A2B3_cI80_199_a2b_2c-001.{AB}                    cI80      I2_{1}3             -         In2O3               #199      In2O3, ICSD #33649
A3B2_cF80_227_f_e-001.{AB}                       cF80      Fd-3m               -         O3Sb2               #227      O3Sb2, ICSD #31102

*****************************************************************************************************************************************************
ANRL - oxides - TERNARIES (R. Friedrich)

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
A4B4C_aP18_2_4i_4i_i-001.{ABC}                   aP18      P-1                 -         Na4O4Si             #2        Na4O4Si, ICSD #15500
A2B7C2_aP22_2_2i_7i_2i-001.{ABC}                 aP22      P-1                 -         Ca2O7V2             #2        Ca2O7V2, ICSD #421266
A2B7C2_aP22_2_2i_7i_2i-002.{ABC}                 aP22      P-1                 -         Mg2O7V2             #2        Mg2O7V2, ICSD #2321
AB3C_aP30_2_3i_9i_3i-001.{ABC}                   aP30      P-1                 -         CaO3Si              #2        CaO3Si, ICSD #201537
A2B5C_aP32_2_4i_10i_2i-001.{ABC}                 aP32      P-1                 -         Al2O5Si             #2        Al2O5Si, ICSD #85742
A2B4C_mP28_4_4a_8a_2a-001.{ABC}                  mP28      P2_{1}              -         Al2O4Sr             #4        Al2O4Sr, ICSD #160297
A3BC_mC60_5_ab8c_ab2c_3c-001.{ABC}               mC60      C2                  -         O3SiSr              #5        O3SiSr, ICSD #32542
A3B5C_mC54_8_3a3b_9a3b_3a-001.{ABC}              mC54      Cm                  -         Ca3O5Si             #8        Ca3O5Si, ICSD #81100
A2B7C3_mP24_11_2e_7e_3e-001.{ABC}                mP24      P2_{1}/m            -         Na2O7Ti3            #11       Na2O7Ti3, ICSD #15463
AB6C2_mC18_12_a_3i_i-001.{ABC}                   mC18      C2/m                -         CaO6V2              #12       CaO6V2, ICSD #166516
AB6C2_mC18_12_a_3i_i-002.{ABC}                   mC18      C2/m                -         MgO6V2              #12       MgO6V2, ICSD #10391
ABC4_mC48_12_gi_hi_2i3j-001.{ABC}                mC48      C2/m                -         MgMoO4              #12       MgMoO4, ICSD #20418
ABC4_mC48_12_gi_hi_2i3j-002.{ABC}                mC48      C2/m                -         FeMoO4;             #12       FeMoO4, ICSD #43012
AB4C_mP12_13_f_2g_e-001.{ABC}                    mP12      P2/c                -         MgO4W               #13       MgO4W, ICSD #67903
AB4C_mP12_13_f_2g_e-002.{ABC}                    mP12      P2/c                -         NiO4W               #13       NiO4W, ICSD #16685
AB3C_mP40_14_2e_6e_2e-001.{ABC}                  mP40      P2_{1}/c            -         MgO3Si              #14       MgO3Si, ICSD #30895
A2B3C_mC24_15_2e_af_e-001.{ABC}                  mC24      C2/c                -         Li2O3Zr             #15       Li2O3Zr, ICSD #94893
AB3C_mC40_15_2e_3f_f-001.{ABC}                   mC40      C2/c                -         NaO3V               #15       NaO3V, ICSD #2103
A2B3C_mC48_15_aef_3f_2e-001.{ABC}                mC48      C2/c                -         Li2O3Ti             #15       Li2O3Ti, ICSD #162215
A4BC7_mC48_15_2f_e_e3f-001.{ABC}                 mC48      C2/c                -         Al4CaO7             #15       Al4CaO7, ICSD #14270
AB3C_mC60_15_cf_e4f_ef-001.{ABC}                 mC60      C2/c                -         CaO3Si              #15       CaO3Si, ICSD #87694
ABC2_oP16_33_a_a_2a-001.{ABC}                    oP16      Pna2_{1}            -         FeNaO2;             #33       FeNaO2, ICSD #186309
A2B3C_oC24_36_b_ab_a-001.{ABC}                   oC24      Cmc2_{1}            -         Li2O3Si             #36       Li2O3Si, ICSD #28192
A2B3C_oC24_36_b_ab_a-002.{ABC}                   oC24      Cmc2_{1}            -         Na2O3Si             #36       Na2O3Si, ICSD #24664
A2B5C_oP32_58_eg_3gh_g-001.{ABC}                 oP32      Pnnm                -         Al2O5Si             #58       Al2O5Si, ICSD #30679
AB2C4_oC28_63_c_ac_fg-001.{ABC}                  oC28      Cmcm                -         CrNa2O4             #63       CrNa2O4, ICSD #76001
AB5C2_oC32_63_c_c2f_f-001.{ABC}                  oC32      Cmcm                -         MgO5Ti2             #63       MgO5Ti2, ICSD #37232
ABC4_tI24_88_a_b_f-001.{ABC}                     tI24      I4_{1}/a            -         BaMoO4              #88       BaMoO4, ICSD #50821
ABC4_tI24_88_a_b_f-002.{ABC}                     tI24      I4_{1}/a            -         CaMoO4              #88       CaMoO4, ICSD #417513
A4BC2_tP28_91_2d_b_ac-001.{ABC}                  tP28      P4_{1}22            -         O4TiZn2             #91       O4TiZn2, ICSD #109093
A2BC4_hP56_173_2b2c_ac_b5c-001.{ABC}             hP56      P6_{3}              -         Al2BaO4             #173      Al2BaO4, ICSD #246027
AB2C4_cF56_227_b_c_e-001.{ABC}                   cF56      Fd-3m               -         MoNa2O4             #227      MoNa2O4, ICSD #151970
    
*****************************************************************************************************************************************************
ANRL - nitrides - BINARIES (R. Friedrich)

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
A2B3_cI80_206_ad_e-001.{ABC}                     cI80      Ia-3                -         N2Zn3               #206      N2Zn3, ICSD #84918

*****************************************************************************************************************************************************
ANRL - nitrides - TERNARIES (R. Friedrich)

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
ABC_oP12_62_c_c_c-001.{ABC}                      oP12      Pnma                -         CaLiN               #62       CaLiN, ICSD #107304
A2B2C_tI10_139_e_e_a-001.{ABC}                   tI10      I4/mmm              -         Ca2N2Zn             #139      Ca2N2Zn, ICSD #69049
A3BC2_hP6_191_f_a_d-001.{ABC}                    hP6       P6/mmm              -         Na3TaTi2            #191      Na3TaTi2, ICSD #186418

*****************************************************************************************************************************************************
ANRL - metal nitrides - BINARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
A2B_mP6_6_3ab_ab-001.{AB}                        mP6       Pm                  -         N2Re1               #6        binary metal-nitride prototype; N2Re1, ICSD #187448
A2B_mC12_8_4a_2a-001.{AB}                        mC12      Cm                  -         N2Re1               #8        binary metal-nitride prototype; N2Re1, ICSD #187447
A3B_mC16_8_6a_b-001.{AB}                         mC16      Cm                  -         N3Rb1               #8        binary metal-nitride prototype; N3Rb1, ICSD #155169
A2B_mP6_11_2e_e-001.{AB}                         mP6       P2_{1}/m            -         N2Re1               #11       binary metal-nitride prototype; N2Re1, ICSD #187444
AB6_mP14_11_e_6e-001.{AB}                        mP14      P2_{1}/m            -         Ba1N6               #11       binary metal-nitride prototype; Ba1N6, ICSD #14244
A3B_mC8_12_di_a-001.{AB}                         mC8       C2/m                -         N3Na1               #12       binary metal-nitride prototype; N3Na1, ICSD #29370
A3B_mC8_12_di_a-002.{AB}                         mC8       C2/m                -         N3Na1               #12       binary metal-nitride prototype; N3Na1, ICSD #29376
AB3_mC8_12_a_di-002.{AB}                         mC8       C2/m                -         Li1N3               #12       binary metal-nitride prototype; Li1N3, ICSD #34675
AB3_mC8_12_a_di-003.{AB}                         mC8       C2/m                -         Li1N3               #12       binary metal-nitride prototype; Li1N3, ICSD #181559
A2B_mC12_12_2i_i-001.{AB}                        mC12      C2/m                -         N2Re1               #12       binary metal-nitride prototype; N2Re1, ICSD #187441
AB_mC16_12_2i_2i-001.{AB}                        mC16      C2/m                -         N1Sr1               #12       binary metal-nitride prototype; N1Sr1, ICSD #411555
A3B2_mC30_12_a4i_3i-001.{AB}                     mC30      C2/m                -         Ca3N2               #12       binary metal-nitride prototype; Ca3N2, ICSD #162794
A3B2_mC30_12_a4i_3i-002.{AB}                     mC30      C2/m                -         Ca3N2               #12       binary metal-nitride prototype; Ca3N2, ICSD #169726
A2B_mP6_13_g_e-001.{AB}                          mP6       P2/c                -         N2Re1               #13       binary metal-nitride prototype; N2Re1, ICSD #187451
AB2_mP12_14_e_2e-001.{AB}                        mP12      P2_{1}/c            -         Ir1N2               #14       binary metal-nitride prototype; Ir1N2, ICSD #160623
A2B_mP12_14_2e_e-004.{AB}                        mP12      P2_{1}/c            -         N2Os1               #14       binary metal-nitride prototype; N2Os1, ICSD #240759
AB3_mP16_14_e_3e-001.{AB}                        mP16      P2_{1}/c            -         Hg2N6               #14       binary metal-nitride prototype; Hg2N6, ICSD #98661
AB2_mC12_15_e_f-001.{AB}                         mC12      C2/c                -         Ba1N2               #15       binary metal-nitride prototype; Ba1N2, ICSD #280681
AB_oC16_36_2a_2a-001.{AB}                        oC16      Cmc2_{1}            -         N1Os1               #36       binary metal-nitride prototype; N1Os1, ICSD #167514
AB_oP8_53_h_h-001.{AB}                           oP8       Pmna                -         N1Os1               #53       binary metal-nitride prototype; N1Os1, ICSD #167513
A2B_oP12_57_e_c-001.{AB}                         oP12      Pbcm                -         N2Re1               #57       binary metal-nitride prototype; N2Re1, ICSD #187450
A2B_oP6_58_g_a-002.{AB}                          oP6       Pnnm                -         N2Os1               #58       binary metal-nitride prototype; N2Os1, ICSD #157283
A2B_oP6_58_g_a-003.{AB}                          oP6       Pnnm                -         N2Pt1               #58       binary metal-nitride prototype; N2Pt1, ICSD #166463
AB2_oP6_58_a_g-004.{AB}                          oP6       Pnnm                -         Ir1N2               #58       binary metal-nitride prototype; Ir1N2, ICSD #160620
AB2_oP6_58_a_g-005.{AB}                          oP6       Pnnm                -         Au1N2               #58       binary metal-nitride prototype; Au1N2, ICSD #166465
A2B_oP12_62_2c_c-005.{AB}                        oP12      Pnma                -         N2Re1               #62       binary metal-nitride prototype; N2Re1, ICSD #187446
AB3_oP16_62_c_cd-002.{AB}                        oP16      Pnma                -         N1Na3               #62       binary metal-nitride prototype; N1Na3, ICSD #165990
AB_oC8_63_c_c-003.AB                             oC8       Cmcm                -         Al1N1               #63       binary metal-nitride prototype; Al1N1, ICSD #163951
A5B3_oC32_63_c2f_cf-001.{AB}                     oC32      Cmcm                -         N5Ta3               #63       binary metal-nitride prototype; N5Ta3, ICSD #16253
A5B3_oC32_63_c2f_cf-002.{AB}                     oC32      Cmcm                -         N5Ta3               #63       binary metal-nitride prototype; N5Ta3, ICSD #66533
A2B_oC12_65_gj_e-001.{AB}                        oC12      Cmmm                -         N2Re1               #65       binary metal-nitride prototype; N2Re1, ICSD #187452
AB_oI4_71_b_c-001.AB                             oI4       Immm                -         Cr1N1               #71       binary metal-nitride prototype; Cr1N1, ICSD #53146
A2B_oI12_71_2e_f-001.{AB}                        oI12      Immm                -         N2Re1               #71       binary metal-nitride prototype; N2Re1, ICSD #187449
AB_oI16_71_abe_n-001.{AB}                        oI16      Immm                -         Li2N2               #71       binary metal-nitride prototype; Li2N2, ICSD #423831
AB3_oI16_72_b_cj-001.{AB}                        oI16      Ibam                -         Ag1N3               #72       binary metal-nitride prototype; Ag1N3, ICSD #27135
A5B4_tI18_87_ah_h-001.{AB}                       tI18      I4/m                -         N5Nb4               #87       binary metal-nitride prototype; N5Nb4, ICSD #26251
AB3_tI32_88_d_cf-001.{AB}                        tI32      I4_{1}/a            -         Cu1N3               #88       binary metal-nitride prototype; Cu1N3, ICSD #30633
AB_tP2_123_d_a-001.AB                            tP2       P4/mmm              -         N1Pr1               #123      binary metal-nitride prototype; N1Pr1, ICSD #168645
A3B_tP4_123_ag_d-001.{AB}                        tP4       P4/mmm              -         N3Rb1               #123      binary metal-nitride prototype; N3Rb1, ICSD #16963
A3B_tP8_123_egh_ab-001.{AB}                      tP8       P4/mmm              -         Cu3N1               #123      binary metal-nitride prototype; Cu3N1, ICSD #180237
A2B_tP6_127_g_b-001.{AB}                         tP6       P4/mbm              -         N2Re1               #127      binary metal-nitride prototype; N2Re1, ICSD #187442
AB_tP4_131_e_c-001.{AB}                          tP4       P4_{2}/mmc          -         N1Tc1               #131      binary metal-nitride prototype; N1Tc1, ICSD #187709
AB_tI4_139_a_b-001.AB                            tI4       I4/mmm              -         Mn1N1               #139      binary metal-nitride prototype; Mn1N1, ICSD #106932
A3B2_tI10_139_ae_e-001.{AB}                      tI10      I4/mmm              -         Mn3N2               #139      binary metal-nitride prototype; Mn3N2, ICSD #71638
A3B2_tI10_139_ae_e-002.{AB}                      tI10      I4/mmm              -         Mn3N2               #139      binary metal-nitride prototype; Mn3N2, ICSD #84202
A3B4_tI14_139_ad_ce-001.{AB}                     tI14      I4/mmm              -         N3Nb4               #139      binary metal-nitride prototype; N3Nb4, ICSD #76389
A8B_tI18_139_deh_a-001.{AB}                      tI18      I4/mmm              -         Fe16N2              #139      binary metal-nitride prototype; Fe16N2, ICSD #41953
AB3_tI16_140_a_dh-001.{AB}                       tI16      I4/mcm              -         K1N3                #140      binary metal-nitride prototype; K1N3, ICSD #24007
AB3_tI16_140_a_dh-002.{AB}                       tI16      I4/mcm              -         Cs1N3               #140      binary metal-nitride prototype; Cs1N3, ICSD #25008
AB3_tI16_140_a_dh-003.{AB}                       tI16      I4/mcm              -         Ag1N3               #140      binary metal-nitride prototype; Ag1N3, ICSD #183201
AB3_tI16_140_a_ch-001.{AB}                       tI16      I4/mcm              -         Cu1N3               #140      binary metal-nitride prototype; Cu1N3, ICSD #187018
AB3_tI16_140_c_ah-001.{AB}                       tI16      I4/mcm              -         K1N3                #140      binary metal-nitride prototype; K1N3, ICSD #1145
AB2_tI12_141_a_e-001.{AB}                        tI12      I4_{1}/amd          -         N1Ti2               #141      binary metal-nitride prototype; N1Ti2, ICSD #23403
A3B_hR4_160_3a_a-001.{AB}                        hR4       R3m                 -         N3Na1               #160      binary metal-nitride prototype; N3Na1, ICSD #644523
A3B2_hR5_160_3a_2a-001.{AB}                      hR5       R3m                 -         Be3N2               #160      binary metal-nitride prototype; Be3N2, ICSD #185490
AB_hP16_162_ek_ci-001.{AB}                       hP16      P-31m               -         Mo1N1               #162      binary metal-nitride prototype; Mo1N1, ICSD #43559
A3B2_hP5_164_ad_d-002.{AB}                       hP5       P-3m1               -         Ca3N2               #164      binary metal-nitride prototype; Ca3N2, ICSD #162795
A3B2_hP5_164_ad_d-003.{AB}                       hP5       P-3m1               -         Ca3N2               #164      binary metal-nitride prototype; Ca3N2, ICSD #162796
A3B2_hP5_164_ad_d-004.{AB}                       hP5       P-3m1               -         Ca3N2               #164      binary metal-nitride prototype; Ca3N2, ICSD #169727
A3B2_hP5_164_ac_d-001.{AB}                       hP5       P-3m1               -         N3V2                #164      binary metal-nitride prototype; N3V2, ICSD #182700
AB_hP16_164_ci_di-001.{AB}                       hP16      P-3m1               -         Mo1N1               #164      binary metal-nitride prototype; Mo1N1, ICSD #60168
A2B_hR3_166_c_a-002.{AB}                         hR3       R-3m                -         Ca2N1               #166      binary metal-nitride prototype; Ca2N1, ICSD #22231
A2B_hR3_166_c_a-003.{AB}                         hR3       R-3m                -         Ba2N1               #166      binary metal-nitride prototype; Ba2N1, ICSD #409851
A3B_hR4_166_bc_a-002.{AB}                        hR4       R-3m                -         N3Na1               #166      binary metal-nitride prototype; N3Na1, ICSD #1144
AB_hP4_186_b_b-002.AB                            hP4       P6_{3}mc            -         Ga1N1               #186      binary metal-nitride prototype; Ga1N1, ICSD #153888
AB_hP4_186_b_b-003.AB                            hP4       P6_{3}mc            -         La1N1               #186      binary metal-nitride prototype; La1N1, ICSD #162195
AB_hP4_186_b_b-004.AB                            hP4       P6_{3}mc            -         Mo1N1               #186      binary metal-nitride prototype; Mo1N1, ICSD #168369
AB_hP4_186_b_a-002.{AB}                          hP4       P6_{3}mc            -         Ga1N1               #186      binary metal-nitride prototype; Ga1N1, ICSD #159250
AB_hP16_186_ac_bc-001.{AB}                       hP16      P6_{3}mc            -         Mo1N1               #186      binary metal-nitride prototype; Mo1N1, ICSD #76280
AB2_hP3_187_a_h-001.{AB}                         hP3       P-6m2               -         Hf1N2               #187      binary metal-nitride prototype; Hf1N2, ICSD #290427
A2B_hP3_191_e_a-001.{AB}                         hP3       P6/mmm              -         N2Os1               #191      binary metal-nitride prototype; N2Os1, ICSD #260545
AB_hP6_191_f_ad-002.{AB}                         hP6       P6/mmm              -         N1Ta1               #191      binary metal-nitride prototype; N1Ta1, ICSD #25659
A2B_hP6_191_af_d-001.{AB}                        hP6       P6/mmm              -         N4Ta2               #191      binary metal-nitride prototype; N4Ta2, ICSD #182351
AB_hP6_191_ad_f-001.{AB}                         hP6       P6/mmm              -         N1Ta1               #191      binary metal-nitride prototype; N1Ta1, ICSD #76457
AB_hP6_191_bc_ad-001.AB                          hP6       P6/mmm              -         N1Y1                #191      binary metal-nitride prototype; N1Y1, ICSD #161079
A3B_hP8_193_g_b-001.{AB}                         hP8       P6_{3}/mcm          -         Ba3N1               #193      binary metal-nitride prototype; Ba3N1, ICSD #77730
AB_hP4_194_c_d-002.AB                            hP4       P6_{3}/mmc          -         Al1N1               #194      binary metal-nitride prototype; Al1N1, ICSD #163950
AB_hP4_194_a_c-001.{AB}                          hP4       P6_{3}/mmc          -         Ir1N1               #194      binary metal-nitride prototype; Ir1N1, ICSD #183157
AB2_hP6_194_c_e-001.{AB}                         hP6       P6_{3}/mmc          -         Hf1N2               #194      binary metal-nitride prototype; Hf1N2, ICSD #290428
A3B2_hP10_194_bf_ac-001.{AB}                     hP10      P6_{3}/mmc          -         Be3N2               #194      binary metal-nitride prototype; Be3N2, ICSD #25656
A3B2_hP10_194_cf_f-001.{AB}                      hP10      P6_{3}/mmc          -         N3W2                #194      binary metal-nitride prototype; N3W2, ICSD #186207
A3B2_hP10_194_bf_ac-002.{AB}                     hP10      P6_{3}/mmc          -         Ca3N2               #194      binary metal-nitride prototype; Ca3N2, ICSD #162797
AB_hP16_194_bh_ag-001.{AB}                       hP16      P6_{3}/mmc          -         Mo1N1               #194      binary metal-nitride prototype; Mo1N1, ICSD #106926
AB3_cI32_204_c_g-001.{AB}                        cI32      Im-3                -         Co1N3               #204      binary metal-nitride prototype; Co1N3, ICSD #162105
AB3_cI32_204_c_g-002.{AB}                        cI32      Im-3                -         Co1N3               #204      binary metal-nitride prototype; Co1N3, ICSD #162106
AB4_cP5_221_a_bc.{AB}                            cP5       Pm-3m               -         N1Ni4               #221      binary metal-nitride prototype; N1Ni4, ICSD #76403

*****************************************************************************************************************************************************
ANRL - metal nitrides - TERNARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
AB2C5_aP8_1_a_2a_5a-001.{ABC}                    aP8       P1                  -         Li1N2Na5            #1        ternary metal-nitride prototype; Li1N2Na5, ICSD #92316
AB2C5_mC16_5_a_c_a2bc-001.{ABC}                  mC16      C2                  -         Li1N2Na5            #5        ternary metal-nitride prototype; Li1N2Na5, ICSD #92315
A3B2C3_mP8_6_2ab_ab_a2b-001.{ABC}                mP8       Pm                  -         Li3N2Na3            #6        ternary metal-nitride prototype; Li3N2Na3, ICSD #92312
A5B2C_mP8_6_3a2b_ab_b-001.{ABC}                  mP8       Pm                  -         Li5N2Na1            #6        ternary metal-nitride prototype; Li5N2Na1, ICSD #92314
A2BC2_mC20_12_2i_i_2i-001.{ABC}                  mC20      C2/m                -         Ca2Fe1N2            #12       ternary metal-nitride prototype; Ca2Fe1N2, ICSD #72389
AB2C12_mC30_12_a_i_6i-001.{ABC}                  mC30      C2/m                -         Cd1K2N12            #12       ternary metal-nitride prototype; Cd1K2N12, ICSD #31297
AB2C12_mC30_12_a_i_6i-002.{ABC}                  mC30      C2/m                -         Cd1K2N12            #12       ternary metal-nitride prototype; Cd1K2N12, ICSD #659621
ABC_mP12_14_e_e_e-001.{ABC}                      mP12      P2_{1}/c            -         Be1Li1N1            #14       ternary metal-nitride prototype; Be1Li1N1, ICSD #402341
ABC_oP12_62_c_c_c-002.{ABC}                      oP12      Pnma                -         Li1Mg1N1            #62       ternary metal-nitride prototype; Li1Mg1N1, ICSD #93259
AB2C_oC16_63_c_2c_c-002.{ABC}                    oC16      Cmcm                -         Au1Ca2N1            #63       ternary metal-nitride prototype; Au1Ca2N1, ICSD #85528
A2BC_oC16_63_2c_c_c-001.{ABC}                    oC16      Cmcm                -         Ca2In1N1            #63       ternary metal-nitride prototype; Ca2In1N1, ICSD #96228
ABC3_oC20_63_c_a_cf-001.{ABC}                    oC20      Cmcm                -         Al1N1Zr3            #63       ternary metal-nitride prototype; Al1N1Zr3, ICSD #29521
A3BC3_oC28_63_cg_c_cg-001.{ABC}                  oC28      Cmcm                -         Ca3Cr1N3            #63       ternary metal-nitride prototype; Ca3Cr1N3, ICSD #40205
A2BC_oC8_65_i_c_b-001.{ABC}                      oC8       Cmmm                -         Li4N2Na2            #65       ternary metal-nitride prototype; Li4N2Na2, ICSD #92307
A7B2C3_oC24_65_a2ij_i_cj-001.{ABC}               oC24      Cmmm                -         Ca7N2Tl3            #65       ternary metal-nitride prototype; Ca7N2Tl3, ICSD #98178
A2BC3_oI12_71_e_a_de-001.{ABC}                   oI12      Immm                -         Ce2Mn1N3            #71       ternary metal-nitride prototype; Ce2Mn1N3, ICSD #50579
ABC_tP6_99_ab_ab_ab-001.{ABC}                    tP6       P4mm                -         Cr1N1Nb1            #99       ternary metal-nitride prototype; Cr1N1Nb1, ICSD #23779
A3B2C3_tP8_115_abc_g_de-001.{ABC}                tP8       P-4m2               -         Li3N2Na3            #115      ternary metal-nitride prototype; Li3N2Na3, ICSD #92311
ABC_tP3_123_c_a_b-001.{ABC}                      tP3       P4/mmm              -         Fe1N1Ni1            #123      ternary metal-nitride prototype; Fe1N1Ni1, ICSD #53505
A5B2C_tP8_123_ai_bc_d-001.{ABC}                  tP8       P4/mmm              -         Li5N2Na1            #123      ternary metal-nitride prototype; Li5N2Na1, ICSD #92313
ABC_tP6_129_c_c_c-001.{ABC}                      tP6       P4/nmm              -         Ca1Ga1N1            #129      ternary metal-nitride prototype; Ca1Ga1N1, ICSD #2027
ABC_tP6_129_c_c_c-002.{ABC}                      tP6       P4/nmm              -         Mo1N1Ta1            #129      ternary metal-nitride prototype; Mo1N1Ta1, ICSD #100437
ABC2_tP8_129_c_c_bc-001.{ABC}                    tP8       P4/nmm              -         Ba1Hf1N2            #129      ternary metal-nitride prototype; Ba1Hf1N2, ICSD #50994
ABC2_tP8_129_a_c_bc-001.{ABC}                    tP8       P4/nmm              -         Li2N2Na4            #129      ternary metal-nitride prototype; Li2N2Na4, ICSD #92309
ABC_tP6_131_e_c_b-001.{ABC}                      tP6       P4_{2}/mmc          -         Ca1N1Ni1            #131      ternary metal-nitride prototype; Ca1N1Ni1, ICSD #69044
ABC_tP6_131_b_c_e-001.{ABC}                      tP6       P4_{2}/mmc          -         Li1N1Sr1            #131      ternary metal-nitride prototype; Li1N1Sr1, ICSD #87414
A2BC_tI8_139_d_b_a-001.ABC/BCA/CAB               tI8       I4/mmm              -         Li4N2Na2            #139      ternary metal-nitride prototype; Li4N2Na2, ICSD #92305
A2B2C_tI10_139_e_e_a-002.{ABC}                   tI10      I4/mmm              -         N2Sr2Zn1            #139      ternary metal-nitride prototype; N2Sr2Zn1, ICSD #80376
A2B2C_tI10_139_e_e_a-003.{ABC}                   tI10      I4/mmm              -         Ba2N2Zn1            #139      ternary metal-nitride prototype; Ba2N2Zn1, ICSD #80377
A3B6C2_tI22_139_ae_eg_e-001.{ABC}                tI22      I4/mmm              -         La3N6V2             #139      ternary metal-nitride prototype; La3N6V2, ICSD #98477
A3B6C2_tI22_139_ae_eg_e-002.{ABC}                tI22      I4/mmm              -         La3N6Nb2            #139      ternary metal-nitride prototype; La3N6Nb2, ICSD #411473
A2B2C_tI20_140_h_h_a-001.{ABC}                   tI20      I4/mcm              -         Be4N4Sr2            #140      ternary metal-nitride prototype; Be4N4Sr2, ICSD #413356
A2BC2_tI20_140_h_a_h-001.{ABC}                   tI20      I4/mcm              -         Be4Ca2N4            #140      ternary metal-nitride prototype; Be4Ca2N4, ICSD #413357
AB2C2_tI20_140_a_h_h-001.{ABC}                   tI20      I4/mcm              -         Ba1Be2N2            #140      ternary metal-nitride prototype; Ba1Be2N2, ICSD #415304
AB2C_hR4_160_a_2a_a-001.{ABC}                    hR4       R3m                 -         Cr1N2W1             #160      ternary metal-nitride prototype; Cr1N2W1, ICSD #84639
A2B2C_hP5_164_d_d_a-001.{ABC}                    hP5       P-3m1               -         Li2N2Zr1            #164      ternary metal-nitride prototype; Li2N2Zr1, ICSD #16231
A2B2C_hP5_164_d_d_a-002.{ABC}                    hP5       P-3m1               -         Mg2N2Sr1            #164      ternary metal-nitride prototype; Mg2N2Sr1, ICSD #410826
AB2C2_hP5_164_a_d_d-001.{ABC}                    hP5       P-3m1               -         Ce1Li2N2            #164      ternary metal-nitride prototype; Ce1Li2N2, ICSD #34003
AB2C_hR4_166_a_c_b-001.{ABC}                     hR4       R-3m                -         Cu1N2Ta1            #166      ternary metal-nitride prototype; Cu1N2Ta1, ICSD #71136
A3BC3_hP14_176_h_c_h-001.{ABC}                   hP14      P6_{3}/m            -         Ba3Fe1N3            #176      ternary metal-nitride prototype; Ba3Fe1N3, ICSD #36502
ABC_hP6_186_b_b_a-001.{ABC}                      hP6       P6_{3}mc            -         N1Na1Sn1            #186      ternary metal-nitride prototype; N1Na1Sn1, ICSD #172471
AB2C_hP8_186_a_ab_b-001.{ABC}                    hP8       P6_{3}mc            -         Mn1N2W1             #186      ternary metal-nitride prototype; Mn1N2W1, ICSD #80029
ABC_hP3_187_a_c_d-001.{ABC}                      hP3       P-6m2               -         Li1N1Ni1            #187      ternary metal-nitride prototype; Li1N1Ni1, ICSD #247028
A5B3C3_hP11_189_dg_g_f-001.{ABC}                 hP11      P-62m               -         Li5N3Ni3            #189      ternary metal-nitride prototype; Li5N3Ni3, ICSD #411152
A2BC_hP4_191_c_a_b-001.{ABC}                     hP4       P6/mmm              -         Li2N1Na1            #191      ternary metal-nitride prototype; Li2N1Na1, ICSD #92308
ABC2_hP4_191_b_a_c-001.{ABC}                     hP4       P6/mmm              -         Li1N1Na2            #191      ternary metal-nitride prototype; Li1N1Na2, ICSD #92310
ABC2_hP8_194_c_a_f-001.{ABC}                     hP8       P6_{3}/mmc          -         Ba1Ce1N2            #194      ternary metal-nitride prototype; Ba1Ce1N2, ICSD #74791
AB2C_hP8_194_a_f_b-001.{ABC}                     hP8       P6_{3}/mmc          -         Fe1N2W1             #194      ternary metal-nitride prototype; Fe1N2W1, ICSD #75971
ABC2_hP8_194_a_c_f-001.{ABC}                     hP8       P6_{3}/mmc          -         Mg1Mo1N2            #194      ternary metal-nitride prototype; Mg1Mo1N2, ICSD #185913
A3BC_hP10_194_h_a_c-001.{ABC}                    hP10      P6_{3}/mmc          -         Ba3N1Na1            #194      ternary metal-nitride prototype; Ba3N1Na1, ICSD #67497
AB3C4_hP16_194_c_ae_2f-001.{ABC}                 hP16      P6_{3}/mmc          -         Al1N3V4             #194      ternary metal-nitride prototype; Al1N3V4, ICSD #181351
AB2C_cF32_227_b_c_a.{ABC}                        cF32      Fd-3m               -         Cs1N2Nb1            #227      ternary metal-nitride prototype; Cs1N2Nb1, ICSD #72546

*****************************************************************************************************************************************************
ANRL - metals - TERNARIES

ANRL Label                                       Pearson   SG symbol           Struk     Prototype           SG #      notes
----------------------------------------------------------------------------------------------------------------------------------------------------
ABC2_mC16_8_b_2a_2ab-001.{ABC}                   mC16      Cm                  -         Bi1Pb1Pd2           #8        metallic ternary prototype; Bi1Pb1Pd2, ICSD #58830
A2B4C3_mP9_10_m_2n_am-001.{ABC}                  mP9       P2/m                -         Co2Gd4Mg3           #10       metallic ternary prototype; Co2Gd4Mg3, ICSD #417035
A2B5C4_mP11_10_m_a2m_2n-001.{ABC}                mP11      P2/m                -         Co2In5Zr4           #10       metallic ternary prototype; Co2In5Zr4, ICSD #55578
AB2C2_mP10_11_e_2e_2e-001.ABC/BCA/CAB            mP10      P2_{1}/m            -         Ba1Bi2Pd2           #11       metallic ternary prototype; Ba1Bi2Pd2, ICSD #416299
A2BC2_mP10_11_2e_e_2e-001.ABC/BCA/CAB            mP10      P2_{1}/m            -         Sn2Sr1Zn2           #11       metallic ternary prototype; Sn2Sr1Zn2, ICSD #424108
A2B4C_mP14_11_2e_4e_e-001.{ABC}                  mP14      P2_{1}/m            -         Au2In4Yb1           #11       metallic ternary prototype; Au2In4Yb1, ICSD #261042
A2BC2_mC10_12_i_a_i-001.{ABC}                    mC10      C2/m                -         Cu2Eu1Sn2           #12       metallic ternary prototype; Cu2Eu1Sn2, ICSD #182050
ABC_mC12_12_i_i_i-001.{ABC}                      mC12      C2/m                -         Al1Ce1Co1           #12       metallic ternary prototype; Al1Ce1Co1, ICSD #20632
A3B2C2_mC14_12_ai_i_i-001.{ABC}                  mC14      C2/m                -         In3Rh2Sr2           #12       metallic ternary prototype; In3Rh2Sr2, ICSD #410985
A6B2C_mC18_12_3i_i_a-001.{ABC}                   mC18      C2/m                -         Na6Sn2Zn1           #12       metallic ternary prototype; Na6Sn2Zn1, ICSD #260159
A2B5C4_mC22_12_i_a2i_2i-001.{ABC}                mC22      C2/m                -         Mn2Sn5Yb4           #12       metallic ternary prototype; Mn2Sn5Yb4, ICSD #419134
A3B4C8_mC30_12_ai_hi_2ij-001.{ABC}               mC30      C2/m                -         Ba3Li4Sn8           #12       metallic ternary prototype; Ba3Li4Sn8, ICSD #240016
A6B2C7_mC30_12_3i_i_a3i-001.{ABC}                mC30      C2/m                -         Ca6Cu2Sn7           #12       metallic ternary prototype; Ca6Cu2Sn7, ICSD #171243
AB2C2_mC20_15_e_f_f-001.{ABC}                    mC20      C2/c                -         In1Pd2Sr2           #15       metallic ternary prototype; In1Pd2Sr2, ICSD #391432
A2B4C_mC28_15_f_2f_e-001.{ABC}                   mC28      C2/c                -         Na2Sn4Sr1           #15       metallic ternary prototype; Na2Sn4Sr1, ICSD #261117
ABC2_oC16_20_b_a_c-001.{ABC}                     oC16      C222_{1}            -         Ag1Bi1K2            #20       metallic ternary prototype; Ag1Bi1K2, ICSD #1156
ABC_oP12_33_a_a_a-001.ABC/BCA/CAB                oP12      Pna2_{1}            -         Ce1Pd1Sn1           #33       metallic ternary prototype; Ce1Pd1Sn1, ICSD #106418
ABC_oC12_36_a_a_a-001.ABC/BCA/CAB                oC12      Cmc2_{1}            -         Ag1Ce1Sn1           #36       metallic ternary prototype; Ag1Ce1Sn1, ICSD #55819
ABC2_oC32_36_b_2a_2ab-001.{ABC}                  oC32      Cmc2_{1}            -         Bi1Pb1Pd2           #36       metallic ternary prototype; Bi1Pb1Pd2, ICSD #56278
ABC_oC18_38_ad_be_bd-001.{ABC}                   oC18      Amm2                -         In1Nd1Rh1           #38       metallic ternary prototype; In1Nd1Rh1, ICSD #59430
ABC_oC18_38_bd_be_ad-001.{ABC}                   oC18      Amm2                -         Cu1Dy1In1           #38       metallic ternary prototype; Cu1Dy1In1, ICSD #260992
ABC_oC18_38_be_ad_bd-001.{ABC}                   oC18      Amm2                -         Gd1In1Pd1           #38       metallic ternary prototype; Gd1In1Pd1, ICSD #415506
A10BC3_oP14_47_2q2rs_h_bt-001.{ABC}              oP14      Pmmm                -         Ga10Ni1Pr3          #47       metallic ternary prototype; Ga10Ni1Pr3, ICSD #20663
ABC2_oP8_51_f_e_be-001.{ABC}                     oP8       Pmma                -         In1La1Ni2           #51       metallic ternary prototype; In1La1Ni2, ICSD #106800
A2BC_oP8_51_ae_f_f-001.{ABC}                     oP8       Pmma                -         Co2Ga1La1           #51       metallic ternary prototype; Co2Ga1La1, ICSD #623100
A2BC_oP8_51_af_f_e-001.{ABC}                     oP8       Pmma                -         Co2Dy1In1           #51       metallic ternary prototype; Co2Dy1In1, ICSD #658675
A4BC_oP12_51_afj_e_e-001.{ABC}                   oP12      Pmma                -         Al4Co1La1           #51       metallic ternary prototype; Al4Co1La1, ICSD #9986
A4BC_oP12_51_afj_e_e-002.{ABC}                   oP12      Pmma                -         Al4Ce1Co1           #51       metallic ternary prototype; Al4Ce1Co1, ICSD #55598
A4BC_oP12_51_afj_e_e-003.{ABC}                   oP12      Pmma                -         In4Rh1Yb1           #51       metallic ternary prototype; In4Rh1Yb1, ICSD #411393
A3B2C_oP12_51_ak_ef_f-001.{ABC}                  oP12      Pmma                -         Au3Rb2Tl1           #51       metallic ternary prototype; Au3Rb2Tl1, ICSD #249924
AB4C_oP12_51_e_afj_e-001.{ABC}                   oP12      Pmma                -         Ca1In4Rh1           #51       metallic ternary prototype; Ca1In4Rh1, ICSD #410891
ABC_oP12_57_c_d_d-001.{ABC}                      oP12      Pbcm                -         Al1Ca1Pd1           #57       metallic ternary prototype; Al1Ca1Pd1, ICSD #370036
AB2C_oP16_57_d_2d_c-001.{ABC}                    oP16      Pbcm                -         Bi1K2Sn1            #57       metallic ternary prototype; Bi1K2Sn1, ICSD #107616
AB2C_oP8_59_a_e_b-001.ABC/BCA/CAB                oP8       Pmmn                -         Cu1Ni2Ti1           #59       metallic ternary prototype; Cu1Ni2Ti1, ICSD #628580
A3B3C_oP14_59_ae_ae_b-001.{ABC}                  oP14      Pmmn                -         Au3In3Sr1           #59       metallic ternary prototype; Au3In3Sr1, ICSD #245679
AB3C4_oP16_59_b_ae_ef-001.{ABC}                  oP16      Pmmn                -         Au1K3Sn4            #59       metallic ternary prototype; Au1K3Sn4, ICSD #107444
AB3C4_oP16_59_b_ae_ef-002.{ABC}                  oP16      Pmmn                -         Au1Cs3Pb4           #59       metallic ternary prototype; Au1Cs3Pb4, ICSD #107448
ABC_oP12_62_c_c_c-003.{ABC}                      oP12      Pnma                -         Ca1Mg1Sn1           #62       metallic ternary prototype; Ca1Mg1Sn1, ICSD #42757
ABC_oP12_62_c_c_c-004.ABC/BCA/CAB                oP12      Pnma                -         Ni1Sn1Tb1           #62       metallic ternary prototype; Ni1Sn1Tb1, ICSD #54301
ABC_oP12_62_c_c_c-005.{ABC}                      oP12      Pnma                -         Co1Dy1Sn1           #62       metallic ternary prototype; Co1Dy1Sn1, ICSD #54415
ABC_oP12_62_c_c_c-006.ABC/BCA/CAB                oP12      Pnma                -         Nd1Pd1Sn1           #62       metallic ternary prototype; Nd1Pd1Sn1, ICSD #54460
ABC_oP12_62_c_c_c-007.{ABC}                      oP12      Pnma                -         Al1Pd1Tb1           #62       metallic ternary prototype; Al1Pd1Tb1, ICSD #54938
ABC_oP12_62_c_c_c-008.{ABC}                      oP12      Pnma                -         Bi1Ca1Li1           #62       metallic ternary prototype; Bi1Ca1Li1, ICSD #58762
ABC_oP12_62_c_c_c-009.{ABC}                      oP12      Pnma                -         Ni1Sc1Sn1           #62       metallic ternary prototype; Ni1Sc1Sn1, ICSD #105338
ABC_oP12_62_c_c_c-010.ABC/BCA/CAB                oP12      Pnma                -         Co1Dy1Sn1           #62       metallic ternary prototype; Co1Dy1Sn1, ICSD #106459
ABC_oP12_62_c_c_c-011.ABC/BCA/CAB                oP12      Pnma                -         La1Ni1Sn1           #62       metallic ternary prototype; La1Ni1Sn1, ICSD #108571
ABC_oP12_62_c_c_c-012.{ABC}                      oP12      Pnma                -         La1Ni1Sn1           #62       metallic ternary prototype; La1Ni1Sn1, ICSD #157921
ABC_oP12_62_c_c_c-013.ABC/BCA/CAB                oP12      Pnma                -         Ni1Pr1Sn1           #62       metallic ternary prototype; Ni1Pr1Sn1, ICSD #157923
ABC_oP12_62_c_c_c-014.{ABC}                      oP12      Pnma                -         Er1Ga1Rh1           #62       metallic ternary prototype; Er1Ga1Rh1, ICSD #630573
AB2C_oP16_62_c_d_c-001.{ABC}                     oP16      Pnma                -         Ca1Cd2Pd1           #62       metallic ternary prototype; Ca1Cd2Pd1, ICSD #425509
AB2C_oP16_62_c_d_c-002.{ABC}                     oP16      Pnma                -         Ga1Pd2Y1            #62       metallic ternary prototype; Ga1Pd2Y1, ICSD #602790
ABC2_oP16_62_c_c_d-001.{ABC}                     oP16      Pnma                -         Dy1Ga1Pd2           #62       metallic ternary prototype; Dy1Ga1Pd2, ICSD #656966
ABC_oC12_63_c_c_c-001.{ABC}                      oC12      Cmcm                -         Mg1Ni1Tb1           #63       metallic ternary prototype; Mg1Ni1Tb1, ICSD #166874
ABC_oC12_63_c_c_a-001.{ABC}                      oC12      Cmcm                -         Pt1Sc1Zn1           #63       metallic ternary prototype; Pt1Sc1Zn1, ICSD #424600
A2BC_oC16_63_f_c_c-001.{ABC}                     oC16      Cmcm                -         Al2Co1Y1            #63       metallic ternary prototype; Al2Co1Y1, ICSD #57645
A2BC_oC16_63_f_c_c-002.{ABC}                     oC16      Cmcm                -         Al2Cu1Mg1           #63       metallic ternary prototype; Al2Cu1Mg1, ICSD #57693
A2BC_oC16_63_f_c_c-003.{ABC}                     oC16      Cmcm                -         Al2Dy1Ni1           #63       metallic ternary prototype; Al2Dy1Ni1, ICSD #57759
A2BC_oC16_63_f_c_c-004.{ABC}                     oC16      Cmcm                -         Ga2Ni1Sc1           #63       metallic ternary prototype; Ga2Ni1Sc1, ICSD #103875
A2BC_oC16_63_f_c_c-005.{ABC}                     oC16      Cmcm                -         In2Rh1Yb1           #63       metallic ternary prototype; In2Rh1Yb1, ICSD #152048
A2BC_oC16_63_f_c_c-006.{ABC}                     oC16      Cmcm                -         Mg2Pd1Sr1           #63       metallic ternary prototype; Mg2Pd1Sr1, ICSD #425493
A2BC_oC16_63_f_c_c-007.{ABC}                     oC16      Cmcm                -         Al2Ce1Pt1           #63       metallic ternary prototype; Al2Ce1Pt1, ICSD #658140
ABC2_oC16_63_a_c_g-001.{ABC}                     oC16      Cmcm                -         Au1Bi1Li2           #63       metallic ternary prototype; Au1Bi1Li2, ICSD #261785
ABC2_oC16_63_a_c_g-002.{ABC}                     oC16      Cmcm                -         Au1Bi1Li2           #63       metallic ternary prototype; Au1Bi1Li2, ICSD #261786
ABC2_oC16_63_a_c_g-003.{ABC}                     oC16      Cmcm                -         Au1Bi1Na2           #63       metallic ternary prototype; Au1Bi1Na2, ICSD #261788
ABC2_oC16_63_a_c_g-004.{ABC}                     oC16      Cmcm                -         Au1Bi1K2            #63       metallic ternary prototype; Au1Bi1K2, ICSD #380341
ABC2_oC16_63_c_c_2c-001.{ABC}                    oC16      Cmcm                -         Ba1Cu1Sn2           #63       metallic ternary prototype; Ba1Cu1Sn2, ICSD #58647
ABC2_oC16_63_c_c_2c-002.{ABC}                    oC16      Cmcm                -         Co1Ho1Sn2           #63       metallic ternary prototype; Co1Ho1Sn2, ICSD #240098
ABC2_oC16_63_c_c_2c-003.{ABC}                    oC16      Cmcm                -         Ce1Ni1Sn2           #63       metallic ternary prototype; Ce1Ni1Sn2, ICSD #621687
ABC2_oC16_63_c_c_2c-004.{ABC}                    oC16      Cmcm                -         Cu1Pr1Sn2           #63       metallic ternary prototype; Cu1Pr1Sn2, ICSD #628742
ABC2_oC16_63_c_c_2c-005.{ABC}                    oC16      Cmcm                -         La1Rh1Sn2           #63       metallic ternary prototype; La1Rh1Sn2, ICSD #641754
AB2C_oC16_63_c_2c_c-003.{ABC}                    oC16      Cmcm                -         Cu1Sn2Sr1           #63       metallic ternary prototype; Cu1Sn2Sr1, ICSD #103110
AB2C_oC16_63_c_2c_c-004.{ABC}                    oC16      Cmcm                -         Co1Sn2Tb1           #63       metallic ternary prototype; Co1Sn2Tb1, ICSD #240096
AB2C_oC16_63_c_f_c-001.{ABC}                     oC16      Cmcm                -         La1Mg2Ni1           #63       metallic ternary prototype; La1Mg2Ni1, ICSD #96152
AB2C_oC16_63_c_f_c-002.{ABC}                     oC16      Cmcm                -         Ce1In2Pd1           #63       metallic ternary prototype; Ce1In2Pd1, ICSD #150132
AB2C_oC16_63_c_f_c-003.{ABC}                     oC16      Cmcm                -         Ba1In2Ir1           #63       metallic ternary prototype; Ba1In2Ir1, ICSD #414434
AB2C_oC16_63_c_f_c-004.{ABC}                     oC16      Cmcm                -         Dy1Ga2Pd1           #63       metallic ternary prototype; Dy1Ga2Pd1, ICSD #600150
A2BC_oC16_63_g_c_c-001.{ABC}                     oC16      Cmcm                -         Cd2Cu1Er1           #63       metallic ternary prototype; Cd2Cu1Er1, ICSD #99139
ABC2_oC16_63_c_c_f-001.{ABC}                     oC16      Cmcm                -         Eu1Pd1Tl2           #63       metallic ternary prototype; Eu1Pd1Tl2, ICSD #172230
A2B2C_oC20_63_g_2c_a-001.{ABC}                   oC20      Cmcm                -         Bi2Cs2Pt1           #63       metallic ternary prototype; Bi2Cs2Pt1, ICSD #658701
A4BC_oC24_63_acf_c_c-001.{ABC}                   oC24      Cmcm                -         Al4Dy1Ni1           #63       metallic ternary prototype; Al4Dy1Ni1, ICSD #57760
A4BC_oC24_63_acf_c_c-002.{ABC}                   oC24      Cmcm                -         Al4Ni1Tb1           #63       metallic ternary prototype; Al4Ni1Tb1, ICSD #58060
AB3C2_oC24_63_c_cg_e-001.{ABC}                   oC24      Cmcm                -         Ca1Ga3Ni2           #63       metallic ternary prototype; Ca1Ga3Ni2, ICSD #58898
A3BC3_oC28_63_cf_a_cf-001.{ABC}                  oC28      Cmcm                -         Co3Ga1Y3            #63       metallic ternary prototype; Co3Ga1Y3, ICSD #10044
ABC2_oC32_63_ac_f_2cf-001.{ABC}                  oC32      Cmcm                -         La1Rh1Sn2           #63       metallic ternary prototype; La1Rh1Sn2, ICSD #410732
A2B2C_oC10_65_j_i_a-001.{ABC}                    oC10      Cmmm                -         Ho2Ni2Pb1           #65       metallic ternary prototype; Ho2Ni2Pb1, ICSD #54612
A2B3C_oC12_65_i_cf_a-001.{ABC}                   oC12      Cmmm                -         Al2Ni3Pr1           #65       metallic ternary prototype; Al2Ni3Pr1, ICSD #107864
A3B2C2_oC14_65_aj_j_i-001.{ABC}                  oC14      Cmmm                -         Mg3Ni2Tb2           #65       metallic ternary prototype; Mg3Ni2Tb2, ICSD #240761
A2BC_oC16_65_aci_j_i-001.{ABC}                   oC16      Cmmm                -         Ga2Nd1Ni1           #65       metallic ternary prototype; Ga2Nd1Ni1, ICSD #103850
A2BC6_oC18_65_j_a_2ij-001.{ABC}                  oC18      Cmmm                -         Dy2Ni1Sn6           #65       metallic ternary prototype; Dy2Ni1Sn6, ICSD #165783
A3B2C7_oC24_65_aj_i_c2ij-001.{ABC}               oC24      Cmmm                -         Ce3Ni2Sn7           #65       metallic ternary prototype; Ce3Ni2Sn7, ICSD #102239
A3B2C7_oC24_65_aj_i_c2ij-002.{ABC}               oC24      Cmmm                -         La3Ni2Sn7           #65       metallic ternary prototype; La3Ni2Sn7, ICSD #160913
A2B2C_oI10_71_f_h_a-001.{ABC}                    oI10      Immm                -         Ce2Ni2Sn1           #71       metallic ternary prototype; Ce2Ni2Sn1, ICSD #55495
A2B2C_oI10_71_f_h_a-002.{ABC}                    oI10      Immm                -         Ca2Cu2Ga1           #71       metallic ternary prototype; Ca2Cu2Ga1, ICSD #58885
A2B2C3_oI14_71_e_f_af-001.{ABC}                  oI14      Immm                -         Al2Sn2Sr3           #71       metallic ternary prototype; Al2Sn2Sr3, ICSD #9564
A5B2C_oI16_71_an_e_d-001.{ABC}                   oI16      Immm                -         Al5Ni2Pr1           #71       metallic ternary prototype; Al5Ni2Pr1, ICSD #58046
A4B3C4_oI22_71_n_af_eg-001.{ABC}                 oI22      Immm                -         Ag8Ce6Sn8           #71       metallic ternary prototype; Ag8Ce6Sn8, ICSD #55826
A4B3C4_oI22_71_n_af_eg-002.{ABC}                 oI22      Immm                -         Ag8Pr6Sn8           #71       metallic ternary prototype; Ag8Pr6Sn8, ICSD #55827
A4B3C4_oI22_71_n_af_eh-001.{ABC}                 oI22      Immm                -         Ag4Dy3Sn4           #71       metallic ternary prototype; Ag4Dy3Sn4, ICSD #156968
A3B4C4_oI22_71_af_eh_n-001.{ABC}                 oI22      Immm                -         La3Pd4Zn4           #71       metallic ternary prototype; La3Pd4Zn4, ICSD #182774
A6B6C_oI26_71_hk_efg_a-001.{ABC}                 oI26      Immm                -         Fe6Ga6Sc1           #71       metallic ternary prototype; Fe6Ga6Sc1, ICSD #103465
A3B3C8_oI28_71_cg_be_2n-001.{ABC}                oI28      Immm                -         Ag3Er3Ga8           #71       metallic ternary prototype; Ag3Er3Ga8, ICSD #605107
A3B3C8_oI28_71_cg_be_2n-002.{ABC}                oI28      Immm                -         Au3Er3Ga8           #71       metallic ternary prototype; Au3Er3Ga8, ICSD #611824
A3B10C_oI28_71_be_g2n_c-001.{ABC}                oI28      Immm                -         Ba3Hg10In1          #71       metallic ternary prototype; Ba3Hg10In1, ICSD #290305
A6B5C3_oI28_71_gn_cn_be-001.{ABC}                oI28      Immm                -         Cu6Sn5Yb3           #71       metallic ternary prototype; Cu6Sn5Yb3, ICSD #413485
A9B3C2_oI28_71_c2n_be_g-001.{ABC}                oI28      Immm                -         Ga9Ho3Pd2           #71       metallic ternary prototype; Ga9Ho3Pd2, ICSD #634397
A8B3C3_oI28_71_2n_be_cg-001.{ABC}                oI28      Immm                -         Ga8Ho3Ir3           #71       metallic ternary prototype; Ga8Ho3Ir3, ICSD #634379
A2B2C_oI20_72_j_j_a-001.{ABC}                    oI20      Ibam                -         Ba2Bi2Zn1           #72       metallic ternary prototype; Ba2Bi2Zn1, ICSD #421424
ABC3_tI10_107_a_a_ab-001.{ABC}                   tI10      I4mm                -         Ba1Ni1Sn3           #107      metallic ternary prototype; Ba1Ni1Sn3, ICSD #58662
A3BC_tI10_107_ab_a_a-001.{ABC}                   tI10      I4mm                -         Al3Cu1Pr1           #107      metallic ternary prototype; Al3Cu1Pr1, ICSD #290390
A3BC_tI10_119_bf_a_c-001.{ABC}                   tI10      I-4m2               -         In3Mg1Sr1           #119      metallic ternary prototype; In3Mg1Sr1, ICSD #249592
AB4C2_tI28_120_d_i_e-001.{ABC}                   tI28      I-4c2               -         Ce1Ni4Sn2           #120      metallic ternary prototype; Ce1Ni4Sn2, ICSD #102238
A4B6C_tI22_121_i_ci_a-001.{ABC}                  tI22      I-42m               -         Ru4Sn6Y1            #121      metallic ternary prototype; Ru4Sn6Y1, ICSD #54354
AB2C_tP4_123_b_h_a-001.ABC/BCA/CAB               tP4       P4/mmm              -         Cd1Pt2Zn1           #123      metallic ternary prototype; Cd1Pt2Zn1, ICSD #102057
AB4C_tP6_123_b_i_a-001.{ABC}                     tP6       P4/mmm              -         Co1Ga4Hf1           #123      metallic ternary prototype; Co1Ga4Hf1, ICSD #623081
ABC5_tP7_123_a_b_ci-001.{ABC}                    tP7       P4/mmm              -         Ce1Co1In5           #123      metallic ternary prototype; Ce1Co1In5, ICSD #102108
AB5C_tP7_123_a_ci_b-001.{ABC}                    tP7       P4/mmm              -         Co1In5Tb1           #123      metallic ternary prototype; Co1In5Tb1, ICSD #623944
AB6C_tP8_123_a_hi_b-001.{ABC}                    tP8       P4/mmm              -         Ce1Ga6Pd1           #123      metallic ternary prototype; Ce1Ga6Pd1, ICSD #240161
AB8C2_tP11_123_a_ehi_g-001.{ABC}                 tP11      P4/mmm              -         Co1Ga8Ho2           #123      metallic ternary prototype; Co1Ga8Ho2, ICSD #42426
AB8C2_tP11_123_a_ehi_g-002.{ABC}                 tP11      P4/mmm              -         Fe1Ga8Ho2           #123      metallic ternary prototype; Fe1Ga8Ho2, ICSD #180132
A2BC2_tP10_127_g_a_h-001.{ABC}                   tP10      P4/mbm              -         Ni2Sn1Zr2           #127      metallic ternary prototype; Ni2Sn1Zr2, ICSD #54303
AB2C2_tP10_127_a_g_h-001.{ABC}                   tP10      P4/mbm              -         Al1Cu2Re2           #127      metallic ternary prototype; Al1Cu2Re2, ICSD #57706
ABC_tP6_129_b_c_c-001.{ABC}                      tP6       P4/nmm              -         Mn1Sn1Sr1           #129      metallic ternary prototype; Mn1Sn1Sr1, ICSD #66951
ABC_tP6_129_c_c_a-001.{ABC}                      tP6       P4/nmm              -         Bi1K1Mn1            #129      metallic ternary prototype; Bi1K1Mn1, ICSD #601586
ABC_tP6_129_a_c_c-001.ABC/BCA/CAB                tP6       P4/nmm              -         Ba1Mg1Pb1           #129      metallic ternary prototype; Ba1Mg1Pb1, ICSD #615958
A2BC_tP8_129_bc_c_a-001.{ABC}                    tP8       P4/nmm              -         Sn2Tb1Zn1           #129      metallic ternary prototype; Sn2Tb1Zn1, ICSD #163425
A2BC_tP8_129_bc_c_a-002.{ABC}                    tP8       P4/nmm              -         Bi2La1Li1           #129      metallic ternary prototype; Bi2La1Li1, ICSD #415728
AB2C2_tP10_129_c_ac_bc-001.ABC/BCA/CAB           tP10      P4/nmm              -         Nd1Ni2Sn2           #129      metallic ternary prototype; Nd1Ni2Sn2, ICSD #160053
AB2C2_tP10_129_c_bc_ac-001.ABC/BCA/CAB           tP10      P4/nmm              -         Eu1Sn2Zn2           #129      metallic ternary prototype; Eu1Sn2Zn2, ICSD #162267
A2BC2_tP10_129_bc_c_ac-001.ABC/BCA/CAB           tP10      P4/nmm              -         Bi2Ce1Ni2           #129      metallic ternary prototype; Bi2Ce1Ni2, ICSD #616562
A2BC4_tP14_129_ac_b_j-001.{ABC}                  tP14      P4/nmm              -         Ce2Ru1Zn4           #129      metallic ternary prototype; Ce2Ru1Zn4, ICSD #418547
AB2C_tI8_139_b_d_a-001.ABC/BCA/CAB               tI8       I4/mmm              -         Li1Pd2Tl1           #139      metallic ternary prototype; Li1Pd2Tl1, ICSD #54365
A2BC2_tI10_139_d_a_e-001.{ABC}                   tI10      I4/mmm              -         Al2Ce1Ga2           #139      metallic ternary prototype; Al2Ce1Ga2, ICSD #55789
A2BC2_tI10_139_d_a_e-002.{ABC}                   tI10      I4/mmm              -         Al2Ca1Zn2           #139      metallic ternary prototype; Al2Ca1Zn2, ICSD #57550
ABC_tI12_139_c_e_e-001.{ABC}                     tI12      I4/mmm              -         Mg1Pr1Sn1           #139      metallic ternary prototype; Mg1Pr1Sn1, ICSD #182479
ABC_tI12_139_e_c_e-001.{ABC}                     tI12      I4/mmm              -         Dy1Mg1Sn1           #139      metallic ternary prototype; Dy1Mg1Sn1, ICSD #183524
ABC_tI12_139_e_c_e-002.ABC/BCA/CAB               tI12      I4/mmm              -         Au1Cu1Sn1           #139      metallic ternary prototype; Au1Cu1Sn1, ICSD #611768
A4B2C_tI14_139_h_d_a-001.{ABC}                   tI14      I4/mmm              -         Al4Mo2Yb1           #139      metallic ternary prototype; Al4Mo2Yb1, ICSD #456
AB4C2_tI14_139_a_h_d-001.{ABC}                   tI14      I4/mmm              -         Er1Ga4Ti2           #139      metallic ternary prototype; Er1Ga4Ti2, ICSD #630579
A5B2C_tI16_139_bg_e_a-001.{ABC}                  tI16      I4/mmm              -         Al5Ni2Zr1           #139      metallic ternary prototype; Al5Ni2Zr1, ICSD #58084
A2BC_tI16_139_ce_e_d-001.{ABC}                   tI16      I4/mmm              -         Bi2Sr1Zn1           #139      metallic ternary prototype; Bi2Sr1Zn1, ICSD #41924
AB2C_tI16_139_e_ce_d-001.{ABC}                   tI16      I4/mmm              -         Ba1Bi2Cd1           #139      metallic ternary prototype; Ba1Bi2Cd1, ICSD #58635
AB2C_tI16_139_e_ce_d-002.{ABC}                   tI16      I4/mmm              -         Ba1Bi2Zn1           #139      metallic ternary prototype; Ba1Bi2Zn1, ICSD #58638
A2BC_tI16_139_ce_d_e-001.{ABC}                   tI16      I4/mmm              -         Bi2Cd1Sr1           #139      metallic ternary prototype; Bi2Cd1Sr1, ICSD #58764
A2B5C3_tI20_139_e_bg_ae-001.{ABC}                tI20      I4/mmm              -         K2Mg5Sn3            #139      metallic ternary prototype; K2Mg5Sn3, ICSD #421342
A8BC4_tI26_139_ij_a_f-001.{ABC}                  tI26      I4/mmm              -         Al8Ca1Cu4           #139      metallic ternary prototype; Al8Ca1Cu4, ICSD #57539
A4B8C_tI26_139_i_fj_a-001.{ABC}                  tI26      I4/mmm              -         Cr4Fe8Y1            #139      metallic ternary prototype; Cr4Fe8Y1, ICSD #168240
A8B4C_tI26_139_ij_f_a-001.{ABC}                  tI26      I4/mmm              -         Al8Mn4Y1            #139      metallic ternary prototype; Al8Mn4Y1, ICSD #57997
ABC4_tI24_140_a_c_l-001.{ABC}                    tI24      I4/mcm              -         Ir1Li1Sn4           #140      metallic ternary prototype; Ir1Li1Sn4, ICSD #172149
A2BC4_tI28_140_h_a_k-001.{ABC}                   tI28      I4/mcm              -         Bi2Mn1Ti4           #140      metallic ternary prototype; Bi2Mn1Ti4, ICSD #150145
A4B2C_tI28_140_l_h_a-001.{ABC}                   tI28      I4/mcm              -         Cu4Sn2Sr1           #140      metallic ternary prototype; Cu4Sn2Sr1, ICSD #182104
AB5C2_tI32_140_a_cl_h-001.{ABC}                  tI32      I4/mcm              -         Bi1Er5Pt2           #140      metallic ternary prototype; Bi1Er5Pt2, ICSD #107217
A3B4C_tI32_140_ah_l_c-001.{ABC}                  tI32      I4/mcm              -         Bi3In4Pb1           #140      metallic ternary prototype; Bi3In4Pb1, ICSD #616725
AB5C2_tI32_140_a_bk_h-001.{ABC}                  tI32      I4/mcm              -         Ga1Nb5Sn2           #140      metallic ternary prototype; Ga1Nb5Sn2, ICSD #103842
A2BC5_tI32_140_h_a_cl-001.{ABC}                  tI32      I4/mcm              -         Au2Bi1Tb5           #140      metallic ternary prototype; Au2Bi1Tb5, ICSD #156956
AB2C5_tI32_140_a_h_cl-001.{ABC}                  tI32      I4/mcm              -         Bi1Co2Ho5           #140      metallic ternary prototype; Bi1Co2Ho5, ICSD #161658
A4B7C3_hR14_146_ab_a2b_b-001.{ABC}               hR14      R3                  -         Ga4Pd7Zn3           #146      metallic ternary prototype; Ga4Pd7Zn3, ICSD #103911
ABC_hP9_156_3a_b2c_2bc-001.{ABC}                 hP9       P3m1                -         Ca1Li1Sn1           #156      metallic ternary prototype; Ca1Li1Sn1, ICSD #58911
A5BC4_hR10_160_5a_a_4a-001.{ABC}                 hR10      R3m                 -         Li5Na1Sn4           #160      metallic ternary prototype; Li5Na1Sn4, ICSD #12142
A2BC3_hP6_164_d_a_bd-002.{ABC}                   hP6       P-3m1               -         Bi2La1Li3           #164      metallic ternary prototype; Bi2La1Li3, ICSD #616769
A2BC_hR4_166_c_b_a-001.{ABC}                     hR4       R-3m                -         Cu2Ga1Sr1           #166      metallic ternary prototype; Cu2Ga1Sr1, ICSD #102938
AB2C_hR4_166_a_c_b-002.{ABC}                     hR4       R-3m                -         Ba1Cu2Ga1           #166      metallic ternary prototype; Ba1Cu2Ga1, ICSD #615828
AB2C4_hR7_166_a_c_2c-001.{ABC}                   hR7       R-3m                -         Au1Ni2Sn4           #166      metallic ternary prototype; Au1Ni2Sn4, ICSD #150127
AB4C6_hR11_166_a_2c_h-001.{ABC}                  hR11      R-3m                -         Gd1Sn4Ti6           #166      metallic ternary prototype; Gd1Sn4Ti6, ICSD #166193
A2BC_hR12_166_h_bc_ac-001.{ABC}                  hR12      R-3m                -         Al2Cu1Yb1           #166      metallic ternary prototype; Al2Cu1Yb1, ICSD #604213
A7B3C2_hR12_166_bh_ac_c-001.{ABC}                hR12      R-3m                -         Al7Ca3Cu2           #166      metallic ternary prototype; Al7Ca3Cu2, ICSD #57538
AB2C9_hR12_166_a_c_bch-001.{ABC}                 hR12      R-3m                -         La1Mg2Ni9           #166      metallic ternary prototype; La1Mg2Ni9, ICSD #55614
A9B2C_hR12_166_eh_c_a-001.{ABC}                  hR12      R-3m                -         Ni3Ti0.67Zr0.33     #166      metallic ternary prototype; Ni3Ti0.67Zr0.33, ICSD #646978
AB2C_hR12_166_bc_h_ac-001.{ABC}                  hR12      R-3m                -         Ag1Al2Pr1           #166      metallic ternary prototype; Ag1Al2Pr1, ICSD #604688
A2B7C3_hR12_166_c_bh_ac-001.{ABC}                hR12      R-3m                -         Ag2Al7Ca3           #166      metallic ternary prototype; Ag2Al7Ca3, ICSD #104173
A7B4C2_hR13_166_ah_2c_c-001.{ABC}                hR13      R-3m                -         Au7Rb4Sn2           #166      metallic ternary prototype; Au7Rb4Sn2, ICSD #58581
A3B2C2_hR14_166_3c_2c_abc-001.{ABC}              hR14      R-3m                -         Au3Sn2Yb2           #166      metallic ternary prototype; Au3Sn2Yb2, ICSD #710044
ABC_hP6_186_b_b_a-002.ABC/BCA/CAB                hP6       P6_{3}mc            -         Cu1Sn1Ti1           #186      metallic ternary prototype; Cu1Sn1Ti1, ICSD #54657
ABC_hP6_186_b_b_a-003.ABC/BCA/CAB                hP6       P6_{3}mc            -         Bi1Cu1Yb1           #186      metallic ternary prototype; Bi1Cu1Yb1, ICSD #416143
ABC_hP6_186_b_a_b-001.ABC/BCA/CAB                hP6       P6_{3}mc            -         Au1Pr1Sn1           #186      metallic ternary prototype; Au1Pr1Sn1, ICSD #54997
ABC_hP6_186_b_a_b-002.{ABC}                      hP6       P6_{3}mc            -         Bi1Li1Zn1           #186      metallic ternary prototype; Bi1Li1Zn1, ICSD #100115
ABC_hP6_186_a_b_b-001.ABC/BCA/CAB                hP6       P6_{3}mc            -         Ca1Sn1Zn1           #186      metallic ternary prototype; Ca1Sn1Zn1, ICSD #106358
A4BC_hP12_186_ac_b_b-001.ABC/BCA/CAB             hP12      P6_{3}mc            -         Cu4In1Mn1           #186      metallic ternary prototype; Cu4In1Mn1, ICSD #424277
A3BC4_hP16_186_a2b_b_2a2b-001.{ABC}              hP16      P6_{3}mc            -         Au3Li1Sn4           #186      metallic ternary prototype; Au3Li1Sn4, ICSD #412207
ABC_hP6_187_i_h_ab-001.ABC/BCA/CAB               hP6       P-6m2               -         Au1Ga1Zr1           #187      metallic ternary prototype; Au1Ga1Zr1, ICSD #156264
A4BC_hP6_187_hi_b_a-001.ABC/BCA/CAB              hP6       P-6m2               -         Ga4Li1Y1            #187      metallic ternary prototype; Ga4Li1Y1, ICSD #98666
ABC_hP9_187_dh_ai_eg-001.{ABC}                   hP9       P-6m2               -         Ca1Li1Pb1           #187      metallic ternary prototype; Ca1Li1Pb1, ICSD #409533
ABC_hP9_187_eh_ci_bg-001.{ABC}                   hP9       P-6m2               -         Li1Pb1Yb1           #187      metallic ternary prototype; Li1Pb1Yb1, ICSD #409534
A5B5C_hP11_187_hk_cgi_a-001.{ABC}                hP11      P-6m2               -         Al5Ba5Pb1           #187      metallic ternary prototype; Al5Ba5Pb1, ICSD #416337
A4B6C_hP11_187_ck_jk_a-001.{ABC}                 hP11      P-6m2               -         Au4In6K1            #187      metallic ternary prototype; Au4In6K1, ICSD #249520
A3B2C3_hP8_189_f_c_g-001.{ABC}                   hP8       P-62m               -         In3Rh2Ti3           #189      metallic ternary prototype; In3Rh2Ti3, ICSD #410967
ABC_hP9_189_g_f_ad-001.{ABC}                     hP9       P-62m               -         Bi1Dy1Rh1           #189      metallic ternary prototype; Bi1Dy1Rh1, ICSD #51845
AB2C6_hP9_189_b_c_fg-001.{ABC}                   hP9       P-62m               -         Co1Ga2Zr6           #189      metallic ternary prototype; Co1Ga2Zr6, ICSD #20876
A5B3C_hP9_189_cf_g_b-001.{ABC}                   hP9       P-62m               -         Ga5Lu3Ni1           #189      metallic ternary prototype; Ga5Lu3Ni1, ICSD #634575
ABC_hP9_189_bc_g_f-001.{ABC}                     hP9       P-62m               -         Ni1Pr1Zn1           #189      metallic ternary prototype; Ni1Pr1Zn1, ICSD #54295
ABC_hP9_189_bc_g_f-002.{ABC}                     hP9       P-62m               -         Pd1Tb1Tl1           #189      metallic ternary prototype; Pd1Tb1Tl1, ICSD #105713
A2BC6_hP9_189_c_b_fg-001.{ABC}                   hP9       P-62m               -         Bi2Fe1Ho6           #189      metallic ternary prototype; Bi2Fe1Ho6, ICSD #96253
ABC_hP9_189_g_bc_f-001.{ABC}                     hP9       P-62m               -         Pb1Pd1Y1            #189      metallic ternary prototype; Pb1Pd1Y1, ICSD #54315
ABC_hP9_189_g_bc_f-002.{ABC}                     hP9       P-62m               -         Hf1Ir1Sn1           #189      metallic ternary prototype; Hf1Ir1Sn1, ICSD #638586
ABC_hP9_189_f_bc_g-001.{ABC}                     hP9       P-62m               -         Mg1Sn1Yb1           #189      metallic ternary prototype; Mg1Sn1Yb1, ICSD #54344
ABC_hP9_189_f_bc_g-002.{ABC}                     hP9       P-62m               -         Ga1Ni1Ti1           #189      metallic ternary prototype; Ga1Ni1Ti1, ICSD #103885
ABC_hP9_189_f_ad_g-001.{ABC}                     hP9       P-62m               -         Dy1In1Mg1           #189      metallic ternary prototype; Dy1In1Mg1, ICSD #55102
ABC_hP9_189_g_f_bc-001.{ABC}                     hP9       P-62m               -         Dy1In1Ni1           #189      metallic ternary prototype; Dy1In1Ni1, ICSD #55846
ABC_hP9_189_f_g_bc-001.{ABC}                     hP9       P-62m               -         Al1Dy1Ni1           #189      metallic ternary prototype; Al1Dy1Ni1, ICSD #107416
ABC_hP9_189_f_g_bc-002.{ABC}                     hP9       P-62m               -         Al1Hf1Pt1           #189      metallic ternary prototype; Al1Hf1Pt1, ICSD #608141
ABC_hP9_189_f_g_bc-003.{ABC}                     hP9       P-62m               -         Al1Nd1Ni1           #189      metallic ternary prototype; Al1Nd1Ni1, ICSD #608770
ABC_hP9_189_ad_g_f-001.{ABC}                     hP9       P-62m               -         Al1Mg1Y1            #189      metallic ternary prototype; Al1Mg1Y1, ICSD #160908
ABC_hP9_189_bc_f_g-001.{ABC}                     hP9       P-62m               -         Cu1In1Tb1           #189      metallic ternary prototype; Cu1In1Tb1, ICSD #628137
A3B2C_hP6_191_g_c_a-001.{ABC}                    hP6       P6/mmm              -         Ag3Al2La1           #191      metallic ternary prototype; Ag3Al2La1, ICSD #57329
A3B2C_hP6_191_g_c_a-002.{ABC}                    hP6       P6/mmm              -         Co3Ga2Tb1           #191      metallic ternary prototype; Co3Ga2Tb1, ICSD #102452
A2BC3_hP6_191_c_a_g-001.{ABC}                    hP6       P6/mmm              -         Cu2Dy1In3           #191      metallic ternary prototype; Cu2Dy1In3, ICSD #108369
A2BC3_hP6_191_c_a_g-002.{ABC}                    hP6       P6/mmm              -         Al2Ce1Pt3           #191      metallic ternary prototype; Al2Ce1Pt3, ICSD #658142
A9BC2_hP12_191_fm_a_c-001.{ABC}                  hP12      P6/mmm              -         Al9Ba1Fe2           #191      metallic ternary prototype; Al9Ba1Fe2, ICSD #57518
A6BC6_hP13_191_i_a_cde-001.{ABC}                 hP13      P6/mmm              -         Mn6Sc1Sn6           #191      metallic ternary prototype; Mn6Sc1Sn6, ICSD #54273
A6B6C_hP13_191_i_cde_a-001.{ABC}                 hP13      P6/mmm              -         Mn6Sn6Tb1           #191      metallic ternary prototype; Mn6Sn6Tb1, ICSD #54277
A3B3C2_hP16_193_g_g_d-001.{ABC}                  hP16      P6_{3}/mcm          -         Ga3Hf3Nb2           #193      metallic ternary prototype; Ga3Hf3Nb2, ICSD #103733
ABC_hP6_194_d_a_c-001.ABC/BCA/CAB                hP6       P6_{3}/mmc          -         Pb1Sr1Zn1           #194      metallic ternary prototype; Pb1Sr1Zn1, ICSD #54319
ABC_hP6_194_a_c_d-001.ABC/BCA/CAB                hP6       P6_{3}/mmc          -         Ag1Ba1Bi1           #194      metallic ternary prototype; Ag1Ba1Bi1, ICSD #56978
ABC_hP6_194_a_c_d-002.ABC/BCA/CAB                hP6       P6_{3}/mmc          -         Ce1Cu1Sn1           #194      metallic ternary prototype; Ce1Cu1Sn1, ICSD #106384
ABC_hP6_194_a_c_d-003.ABC/BCA/CAB                hP6       P6_{3}/mmc          -         Ce1Tl1Zn1           #194      metallic ternary prototype; Ce1Tl1Zn1, ICSD #622301
ABC_hP6_194_c_a_d-001.ABC/BCA/CAB                hP6       P6_{3}/mmc          -         Bi1Ca1Cu1           #194      metallic ternary prototype; Bi1Ca1Cu1, ICSD #57018
ABC_hP6_194_c_a_d-002.ABC/BCA/CAB                hP6       P6_{3}/mmc          -         Cu1Gd1Sn1           #194      metallic ternary prototype; Cu1Gd1Sn1, ICSD #600994
ABC_hP6_194_a_d_c-001.ABC/BCA/CAB                hP6       P6_{3}/mmc          -         Ba1Bi1Cu1           #194      metallic ternary prototype; Ba1Bi1Cu1, ICSD #106303
ABC_hP6_194_a_d_c-002.ABC/BCA/CAB                hP6       P6_{3}/mmc          -         Ca1Hg1Pb1           #194      metallic ternary prototype; Ca1Hg1Pb1, ICSD #106345
ABC_hP6_194_c_d_a-002.ABC/BCA/CAB                hP6       P6_{3}/mmc          -         Hg1Pb1Yb1           #194      metallic ternary prototype; Hg1Pb1Yb1, ICSD #108545
AB2C_hP8_194_b_f_c-001.{ABC}                     hP8       P6_{3}/mmc          -         Al1Pt2Zr1           #194      metallic ternary prototype; Al1Pt2Zr1, ICSD #58139
AB2C_hP8_194_a_f_c-001.{ABC}                     hP8       P6_{3}/mmc          -         In1Pt2Y1            #194      metallic ternary prototype; In1Pt2Y1, ICSD #59505
AB2C_hP8_194_c_f_a-001.{ABC}                     hP8       P6_{3}/mmc          -         Er1Pt2Sn1           #194      metallic ternary prototype; Er1Pt2Sn1, ICSD #103292
AB2C_hP8_194_c_f_a-002.{ABC}                     hP8       P6_{3}/mmc          -         Ce1Cu2Mg1           #194      metallic ternary prototype; Ce1Cu2Mg1, ICSD #157110
A2BC_hP8_194_f_a_c-001.{ABC}                     hP8       P6_{3}/mmc          -         Cu2Li1Sn1           #194      metallic ternary prototype; Cu2Li1Sn1, ICSD #150602
A3BC2_hP12_194_h_a_f-001.{ABC}                   hP12      P6_{3}/mmc          -         Al3Ru1Sc2           #194      metallic ternary prototype; Al3Ru1Sc2, ICSD #58159
ABC_hP12_194_f_f_ab-001.{ABC}                    hP12      P6_{3}/mmc          -         Ga1Sn1Sr1           #194      metallic ternary prototype; Ga1Sn1Sr1, ICSD #66003
ABC_hP12_194_ab_f_f-001.{ABC}                    hP12      P6_{3}/mmc          -         Gd1Sn1Zn1           #194      metallic ternary prototype; Gd1Sn1Zn1, ICSD #152621
ABC_hP12_194_ab_f_f-002.{ABC}                    hP12      P6_{3}/mmc          -         La1Sn1Zn1           #194      metallic ternary prototype; La1Sn1Zn1, ICSD #152623
ABC_hP12_194_f_ab_f-001.{ABC}                    hP12      P6_{3}/mmc          -         Sn1Y1Zn1            #194      metallic ternary prototype; Sn1Y1Zn1, ICSD #152622
A3B3C2_hP16_194_h_af_f-001.{ABC}                 hP16      P6_{3}/mmc          -         Al2Ba2Ga1.33        #194      metallic ternary prototype; Al2Ba2Ga1.33, ICSD #415931
AB2C_cP16_215_e_bce_ad-001.{ABC}                 cP16      P-43m               -         Ir1Li2Mg1           #215      metallic ternary prototype; Ir1Li2Mg1, ICSD #180119
A3B4C4_cI22_217_b_c_c-001.{ABC}                  cI22      I-43m               -         Ca3Ga4Ni4           #217      metallic ternary prototype; Ca3Ga4Ni4, ICSD #58899
A3BC_cP15_221_ag_c_d-001.{ABC}                   cP15      Pm-3m               -         Ga3Os1Tb1           #221      metallic ternary prototype; Ga3Os1Tb1, ICSD #412144
AB3C_cP15_221_d_ag_c-001.{ABC}                   cP15      Pm-3m               -         Dy1Ga3Ru1           #221      metallic ternary prototype; Dy1Ga3Ru1, ICSD #629734

****************************************************************************************************
****************************************************************************************************
****************************************************************************************************
****************************************************************************************************
*                                                                                                  *
*                        aflow - STEFANO CURTAROLO Duke University 2003-2020                       *
*                            High-Throughput ab-initio Computing Project                           *
*                                                                                                  *
****************************************************************************************************
