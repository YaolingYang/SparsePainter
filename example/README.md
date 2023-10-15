This example includes 8000 reference individuals from 4 populations with 2091 SNPs (Both vcf version ``donor.vcf.gz`` and phase version ``donor.phase.gz`` are available), and the aim is to paint 500 target individuals (Both vcf version ``target.vcf.gz`` and phase version ``target.phase.gz`` are available). Remember we have compiled SparsePainter in ``SparsePainter``, then we can paint with the following command:

(a) If your input file is in vcf or vcf.gz format:

``
./SparsePainter -reffile donor.vcf.gz -targetfile target.vcf.gz -popfile popnames.txt -mapfile map.txt -namefile targetname.txt -out target_vs_ref -prob -chunklength -aveSNP -aveind
``

(b) If your input file is in phase of phase.gz format:

``
./SparsePainter -reffile donor.phase.gz -targetfile target.phase.gz -popfile popnames.txt -mapfile map.txt -namefile targetname.txt -out target_vs_ref -prob -chunklength -aveSNP -aveind
``

The output file for this example includes ``target_vs_ref_prob.txt.gz``, ``target_vs_ref_chunklength.txt.gz``, ``target_vs_ref_aveSNPprob.txt``, ``target_vs_ref_aveindprob.txt`` and ``target_vs_ref_fixedlambda.txt``.

To paint the reference individuals against themselves with leave-one-out strategy, run with:

(a) If your input file is in vcf or vcf.gz format:

``
./SparsePainter -reffile donor.vcf.gz -targetfile donor.vcf.gz -popfile popnames.txt -mapfile map.txt -namefile refname.txt -out ref_vs_ref -prob -chunklength -aveSNP -aveind -loo
``

(b) If your input file is in phase or phase.gz format:

``
./SparsePainter -reffile donor.phase.gz -targetfile donor.phase.gz -popfile popnames.txt -mapfile map.txt -namefile refname.txt -out ref_vs_ref -prob -chunklength -aveSNP -aveind -loo
``

The output file for this example includes ``ref_vs_ref_prob.txt.gz``, ``ref_vs_ref_chunklength.txt.gz``, ``ref_vs_ref_aveSNPprob.txt``, ``ref_vs_ref_aveindprob.txt`` and ``ref_vs_ref_fixedlambda.txt``.
