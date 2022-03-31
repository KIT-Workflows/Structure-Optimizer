import os,re, tarfile, yaml, tarfile


if __name__ == '__main__':

    
    conf_algo = {'iMTD-GC':"-v3",'MF-MD-GC':"-v1", 'iMTD-sMTD':"-v4", 'Conformational entropy':"-entropy"}

    with open('rendered_wano.yml') as file:
        wano_file = yaml.full_load(file)

    solv_model = wano_file["Parameters"]["Solvent"]
    conf_method = conf_algo[ wano_file["Parameters"]["Conformational-algo"] ]
    Opt_level = wano_file["Parameters"]["Optimization level"]
    e_thr = str(float(wano_file["Parameters"]["E-thr"]))
    r_thr = str(float(wano_file["RMSD-thr"]))


    with open('Conformer-gen.slr') as file:
        lines = file.readlines()
    
    n_threads = re.findall(r"\w+", lines[4])[4]

    input_file = wano_file["Input-File"]

    path_CREST = '/home/ws/gt5111/XTB/crest-2.11.2/_build_intel/crest '

    # Create bash file
    bash_file_name = "run_crest.sh"
    with open(bash_file_name, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('. /etc/profile.d/lmod.sh\n')
        f.write('set -e\n')

        f.write('module purge\n')
        f.write('module load xtb/6.4.1\n')
        f.write('module load intel/19.0.5.281\n')

        f.write('#defaults\n')
        f.write('solv="gas"\n')   # solvent name
        f.write('ewin="8.0"\n')   # energy window, do not take < 8 for larger drugs
        f.write('gbsa=' + "-" + solv_model + '\n') # solvation model
        f.write('mode="ff"\n')    # FF+GFN2(-screen)

        # settings for full, non-quick default mode, change by -q 
        
        f.write('listd="0.0 0.5 1.5 2.0"\n') # artificial dispersion
        f.write('listc="1 2 -1 -2"\n')       #    "       ES
        f.write('xtff="x1.0"\n')             # neutral run time
        f.write('xtffm="x2.0"\n')            # artifical "  "
        f.write('nclust="100"\n')            # of clusters in itermediate steps
        f.write('vfirst=' + conf_method + '\n')            # algo for first neutral run
        f.write('nmr=""\n')                  # create anmr_nucinfo and anmr_rotamer in first FF(v4) run

        #f.write('threads=$OMP_NUM_THREADS\n') # set OMP threads (equivalent to using the -T option of crest)
        f.write('threads=' + n_threads + '\n')
        f.write('printhelp="off"\n')
        f.write('argcheck="y"\n')       

        f.write('while [ "$argcheck" = "y" ]; do\n')
        f.write('   if [ -n "$1" ]; then\n')
        f.write('      case $1 in\n')
        f.write('      "-q"    )           listd="0.5 1.5"; listc="1 -1"; xtffm="x0.5"; vfirst="-v3"; mode="ff"; nclust="400" ;;\n')
        f.write('      "-s"    ) shift;    solv=$1 ;;\n')
        f.write('      "-solvent" ) shift; solv=$1 ;;\n')
        f.write('      "-g"    ) shift;    solv=$1 ;;\n')
        f.write('      "-ewin" ) shift;    ewin=$1 ;;\n')
        f.write('      "-xt" )             mode="" ;;\n')    #xt = extended search with FF + GFN2 PES search
        f.write('      "-T" ) shift;       threads=$1 ;;\n')
        f.write('      "-threads" ) shift; threads=$1 ;;\n')
        f.write('      "-nmr" )            nmr="-nmr" ;;\n')
        f.write('      "-or")              nclust="1000";;\n')
        f.write('      "-h"      ) printhelp="on";;\n')
        f.write('      "--help"  ) printhelp="on";;\n')
        f.write('      esac\n')
        f.write('      shift\n')
        f.write('   else\n')
        f.write('      argcheck="n"\n')
        f.write('   fi\n')  
        f.write('  done\n')        

        f.write('ewin2=$(echo "$ewin + 4" |bc)\n') # ewin for intermediate FF run re-ranking with GFN2
        f.write('ewin3=$(echo "$ewin - 2" |bc)\n') # ewin for intermediate FF run re-ranking with GFNlast GFN2 search (if mode="")
        
        f.write('echo ""\n')
        f.write('echo "*****************************************"\n')
        f.write('echo "CREST energy window      : $ewin kcal/mol"\n')
        f.write('echo "ewin for FF reranking    : $ewin2 kcal/mol"\n')
        f.write('echo "applied solvent          : $solv"\n')
        f.write('if [ "$mode" == "" ]; then\n')
        f.write('    echo "search method            : FF + GFN2 "\n')
        f.write('    echo "ewin for last GFN2 search: $ewin3"\n')
        f.write('else\n')
        f.write('    echo "search method            : $mode"\n')
        f.write('fi\n')
        f.write('echo "nclust                   : $nclust"\n')
        f.write('echo "*****************************************"\n')

        #########################
        #search on FF level
        #########################
        # neutral
        f.write('rm -f gfnff_topo\n')
        f.write('echo -e "#>>># Running: CREST 1 (FF,$vfirst)" >&2\n')
        f.write('echo -e "#>>># Running: CREST 1 (FF,$vfirst)"\n') 

        f.write('if [ "$nmr" == "-nmr" ]; then\n')
        f.write('    echo "Generating anmr_nucinfo and anmr_rotamers files for possible subsequent NMR calculation"\n')
        f.write(path_CREST + ' ' + input_file + ' "$gbsa" "$solv" -gfnff "$vfirst" -cluster "$nclust" -ewin "$ewin" -norotmd -nocross -mdlen "$xtff" -T "$threads" "$nmr"\n')
        f.write('else\n')
        f.write(path_CREST + ' '+ input_file + ' "$gbsa" "$solv" -gfnff "$vfirst" -cluster "$nclust" -ewin "$ewin" -norotmd -nocross -mdlen "$xtff" -T "$threads"\n')
        f.write('fi\n')
        # move old files
        f.write('if [ -f "crest_conformers.xyz" ]; then\n')
        f.write('    mv crest_conformers.xyz crest_tmp.xyz\n')
        f.write('fi\n')
        f.write('if [ -f "crest_best.xyz" ]; then\n')
        f.write('    mv crest_best.xyz best_gfnff.xyz\n')
        f.write('fi\n')
        f.write('if [ -f "crest_clustered.xyz" ]; then\n')
        f.write('    cat crest_clustered.xyz > search_gfnff.xyz\n')
        f.write('fi\n')

        # artificial 1
        f.write('for dscal in $listd\n')
        f.write('do\n')
        f.write('echo -e "#>>># Running: CREST (FF+/-disp)... $dscal" >&2\n')
        f.write('echo -e "#>>># Running: CREST (FF+/-disp)... $dscal"\n')
        f.write('rm -f gfnff_topo\n')
        # search with scaled dispersion interaction via the "-dispscal" command
        f.write(path_CREST + ' best_gfnff.xyz -gfnff -cluster "$nclust" -ewin "$ewin" -dispscal "$dscal" -norotmd -norestart -nocross -mdlen "$xtffm" -T "$threads"\n')
        f.write('cat crest_clustered.xyz >> search_gfnff.xyz\n')
        f.write('done\n')

        f.write('chrg=0\n')
        f.write('if test -f .CHRG; then\n')
        f.write('chrg=$(cat .CHRG)\n')
        f.write('fi\n')

        # artificial 2
        f.write('for c in $listc\n')
        f.write('do\n')
        f.write('chrgmod=$((chrg+c))\n')
        f.write('echo -e "#>>># Running: CREST (FF+/-charge)... $chrgmod" >&2\n')
        f.write('echo -e "#>>># Running: CREST (FF+/-charge)... $chrgmod"\n')
        #echo 'CREST (FF+/-charge)...' $chrgmod
        f.write('rm -f gfnff_topo\n')
        f.write('echo "$chrgmod" > .CHRG\n')
        # search with modified molecular charge +/-{1,2}
        f.write(path_CREST + ' best_gfnff.xyz -gfnff -cluster "$nclust" -ewin "$ewin" -norotmd -norestart -nocross -mdlen "$xtffm" -T "$threads"\n')
        f.write('cat crest_clustered.xyz >> search_gfnff.xyz\n')
        f.write('done\n')
        f.write('echo "$chrg" > .CHRG\n')

        # GFN2 opt
        f.write('echo "CREST screen (GFN2)..."\n') 
        f.write('cp search_gfnff.xyz search_gfnff_save.xyz\n')
        f.write('rm -f gfnff_topo\n')
        #1. energy re-rank (at FF level) because of +/-disp/ES structures (output: crest_ensemble.xyz)
        f.write('echo -e "#>>># Running: CREST structure optimization and energy reranking because of artificial search" >&2\n')
        f.write('echo -e "#>>># Running: CREST structure optimization and energy reranking because of artificial search"\n')
        f.write(path_CREST + ' -mdopt search_gfnff.xyz "$gbsa" "$solv" -opt crude -gfnff -ewin "$ewin2" -T "$threads"\n')  
        #2. Sorting of the optimized ensemble (output: overwritten crest_ensemble.xyz)
        f.write('echo -e "#>>># Running: CREST (cregen) ensemble sorting" >&2\n')
        f.write('echo -e "#>>># Running: CREST (cregen) ensemble sorting"\n')
        f.write(path_CREST + ' -cregen crest_ensemble.xyz -ewin "$ewin" -ethr ' + e_thr + ' und -rthr ' + r_thr +' -bthr 0.03 -T "$threads"\n')
        #3. GFN2 full opt (output: overwritten crest_ensemble.xyz)
        f.write('echo -e "#>>># Running: CREST GFN2 full optimization" >&2\n')
        f.write('echo -e "#>>># Running: CREST GFN2 full optimization"\n')
        f.write(path_CREST + ' "$gbsa" "$solv" -screen crest_ensemble.xyz -opt ' + Opt_level + ' -gfn2 -ewin "$ewin" -T "$threads"\n')
        #4. PCA/k-Means clustering of final ensemble (output: crest_clustered.xyz)
        f.write('echo -e "#>>># Running: CREST initial FF PCA/k-Means clustering" >&2\n')
        f.write('echo -e "#>>># Running: CREST initial FF PCA/k-Means clustering"\n')
        f.write(path_CREST + ' -for crest_ensemble.xyz -cluster "$nclust" -ewin "$ewin" -T "$threads"\n')
        f.write('mv crest_clustered.xyz search_gfnff.xyz\n')
        f.write('mv crest_best.xyz        best_gfnff.xyz\n')

        #########################
        # FF only part ends here
        #########################
        f.write('if [ "$mode" = "ff" ]; then\n')
        f.write('    echo -e "#>>># Running: CREST final PCA/k-Means clustering for FF" >&2\n')
        f.write('    echo -e "#>>># Running: CREST final PCA/k-Means clustering for FF"\n') 
        f.write(path_CREST + ' -for search_gfnff.xyz -cluster tightincr -ewin "$ewin" -T "$threads"\n') 
        f.write('    mv crest_clustered.xyz crest_combi.xyz\n')
        f.write('else\n')

        ###############################################################
        #search on GFN2 level starting with best at GFN2 opt from FF SE
        ###############################################################
        f.write('echo "CREST 8 (GFN2) ... "\n')
        f.write('echo -e "#>>># Running: CREST GFN2 search starting from FF SE = CREST 8 (GFN2) " >&2\n')
        f.write('echo -e "#>>># Running: CREST GFN2 search starting from FF SE = CREST 8 (GFN2) "\n')
        f.write('crest best_gfnff.xyz "$gbsa" "$solv" -gfn2 -ewin "$ewin3" -norotmd -nocross -mdlen x1.0 -T "$threads"\n')
        f.write('echo -e "#>>># Running: CREST GFN2 PCA/k-Means clustering into $nclust cluster" >&2\n')
        f.write('echo -e "#>>># Running: CREST GFN2 PCA/k-Means clustering into $nclust cluster"\n')
        f.write(path_CREST + ' -for crest_conformers.xyz -cluster "$nclust" -T "$threads"\n')
        f.write('mv crest_best.xyz        best_gfn2.xyz\n')
        f.write('mv crest_clustered.xyz search_gfn2.xyz\n')
        
        # merge FF/TB  
        f.write('cat search_gfnff.xyz   > final.xyz\n')
        f.write('cat search_gfn2.xyz   >> final.xyz\n')
        f.write('echo "CREST final cregen ..."\n')
        f.write('echo -e "#>>># Running: CREST (cregen) sorting final (FF/SQM) ensemble" >&2\n')
        f.write('echo -e "#>>># Running: CREST (cregen) sorting final (FF/SQM) ensemble"\n')
        f.write(path_CREST + ' -cregen final.xyz -ewin "$ewin" -T "$threads"\n')
        f.write('echo -e "#>>># Running: CREST PCA/k-Means clustering for final (FF/SQM) ensemble" >&2\n')
        f.write('echo -e "#>>># Running: CREST PCA/k-Means clustering for final (FF/SQM) ensemble"\n')
        f.write(path_CREST + ' -for crest_ensemble.xyz -cluster "$nclust" -ewin "$ewin" -T "$threads"\n')
        f.write('echo -e "#>>># Running: CREST final PCA/k-Means clustering for final (FF/SQM) ensemble" >&2\n')
        f.write('echo -e "#>>># Running: CREST final PCA/k-Means clustering for final (FF/SQM) ensemble"\n')
        f.write(path_CREST + ' -for crest_clustered.xyz -cluster tightincr -ewin "$ewin" -T "$threads"\n')
        f.write('cp crest_clustered.xyz crest_combi.xyz\n')
        f.write('fi\n')

        # the workflow provides a final ensemble called "crest_combi.xyz" for further processing

        f.write('echo -e "#>>># CREST_COMBI: All done!" >&2\n')
        f.write('echo "#>>># CREST_COMBI: All done!"\n')
        f.write('\n')

    os.system("chmod +x " + bash_file_name)


    