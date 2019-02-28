import os
import os.path
import shutil
import subprocess
import sys

PARAM_BINARY = './JTC/cc/parameter'
DATA_DIR = 'tmp/real_input_data/'
OUTPUT_DIR = 'tmp/real_data_runs/'

def runExperiment(args):
        original_data_file_path = DATA_DIR + args['data_file']
        if not os.path.isfile(original_data_file_path):
            sys.stderr.write("Data file doesn't exist: %s\n" % original_data_file_path)
            sys.exit(2)
        output_dir = os.path.join(OUTPUT_DIR, str(args['run']))

        os.makedirs(output_dir)
        data_file_path = os.path.join(output_dir, "data")
        shutil.copyfile(original_data_file_path, data_file_path)

        output_file_path = os.path.join(output_dir, 'output')
        output_file = file(output_file_path, 'w')
        call_arguments = [PARAM_BINARY, '--noapproximate', '-f', data_file_path, '-i', str(args['iterations'])]
        if 'ref_csm' in args:
            call_arguments.append("--ref_csm")
            call_arguments.append(str(args['ref_csm']))
        if 'ref_csv' in args:
            call_arguments.append("--ref_csv")
            call_arguments.append(str(args['ref_csv']))
        if 'ref_tm' in args:
            call_arguments.append("--ref_tm")
            call_arguments.append(str(args['ref_tm']))
        if 'ref_tv' in args:
            call_arguments.append("--ref_tv")
            call_arguments.append(str(args['ref_tv']))
        print "Calling prog: " + ' '.join(call_arguments)
        try:
            subprocess.call(
                    call_arguments,
                    stdout=output_file)
        finally: 
            output_file.close()

if __name__ == "__main__":
    params = [
        {'run': 0, 'iterations':25,
                'data_file': 'block1/all_data',
                'ref_csm': -2.0, 'ref_csv': 17.0,
                'ref_tm': 3.0, 'ref_tv': 5},
        {'run': 1, 'iterations':25,
                'data_file': 'block1/control_data',
                'ref_csm': -2.0, 'ref_csv': 17.0,
                'ref_tm': 3.0, 'ref_tv': 5},
        {'run': 2, 'iterations':25,
                'data_file': 'block1/unhealthy_data',
                'ref_csm': -2.0, 'ref_csv': 17.0,
                'ref_tm': 3.0, 'ref_tv': 5},
        {'run': 3, 'iterations':25,
                'data_file': 'block2/all_data',
                'ref_csm': -2.0, 'ref_csv': 17.0,
                'ref_tm': 3.0, 'ref_tv': 5},
        {'run': 4, 'iterations':25,
                'data_file': 'block2/control_data',
                'ref_csm': -2.0, 'ref_csv': 17.0,
                'ref_tm': 3.0, 'ref_tv': 5},
        {'run': 5, 'iterations':25,
                'data_file': 'block2/unhealthy_data',
                'ref_csm': -2.0, 'ref_csv': 17.0,
                'ref_tm': 3.0, 'ref_tv': 5},
        {'run': 6, 'iterations':25,
                'data_file': 'block3/all_data',
                'ref_csm': -2.0, 'ref_csv': 17.0,
                'ref_tm': 3.0, 'ref_tv': 5},
        {'run': 7, 'iterations':25,
                'data_file': 'block3/control_data',
                'ref_csm': -2.0, 'ref_csv': 17.0,
                'ref_tm': 3.0, 'ref_tv': 5},
        {'run': 8, 'iterations':25,
                'data_file': 'block3/unhealthy_data',
                'ref_csm': -2.0, 'ref_csv': 17.0,
                'ref_tm': 3.0, 'ref_tv': 5},
    ]
    if len(sys.argv) > 2:
        sys.stderr.write("Invalid arguments. Usage: python real_data_experiment.py [experiment_id]\n")
        sys.exit(1)
    elif len(sys.argv) == 2:
        experiment_id = int(sys.argv[1])
        runExperiment(params[experiment_id])
    else:
        for args in params:
            runExperiment(args)
