import os, tarfile, shutil, yaml

def gen_xyz_files(file_name):
    
    with open(file_name) as f:
        first_line = f.readline()

    lines_per_file = int(first_line) + int(2)
    smallfile = None

    with open(file_name) as bigfile:
        for lineno, line in enumerate(bigfile):
            if lineno % lines_per_file == 0:
                if smallfile:
                    smallfile.close()
                small_filename = 'mol_{}.xyz'.format(lineno + lines_per_file)
                smallfile = open(small_filename, "w")
            smallfile.write(line)
        if smallfile:
            smallfile.close()
    
    return print("The .xyz files were successfully written!")

if __name__ == '__main__':

    with open('rendered_wano.yml') as file:
        wano_file = yaml.full_load(file)

    engine = wano_file = wano_file["Parameters"]["Engine"]

    if engine == "RdKit":
        os.system('python conf_rdkit.py')
    else:
        os.system('python crest.py')
        os.system('bash run_crest.sh')
        
        file_name = 'crest_combi.xyz'
        archive = tarfile.open("mol.tar.gz", "w|gz")

        gen_xyz_files(file_name)

        for fname in os.listdir():
            if fname.startswith("mol_"):
                archive.add(fname, arcname=fname)
                os.remove(os.path.join(fname)) 
        archive.close()

        # file_name = 'crest_combi.xyz'

    # archive = tarfile.open("mol.tar.gz", "w|gz")

    # gen_xyz_files(file_name)

    # for fname in os.listdir():
    #     if fname.startswith("mol_"):
    #         archive.add(fname, arcname=fname)
    #         os.remove(os.path.join(fname)) 

    # archive.close()                    

    # input_file = wano_file["Input-File"]
    # mutability = wano_file["mutability"]
    # nconf = wano_file["nconf"]
    

    # os.system('obabel Henrik12.xyz -O conformers.xyz --confab --conf 50')
    # os.system('obabel conformers.xyz -oconfabreport -xf Henrik12.xyz -xr 1.0')





