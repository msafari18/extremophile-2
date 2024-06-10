import argparse
import os
from deep_clustering import run
from build_signature_dataset import run_fragment_builder



ENVS = ["Radio_temp"]
NUM_CLUSTERS = {"Radio_temp": [2, 4]}
# FRAGMENT_LENGTHS = [10000, 50000, 100000, 250000, 500000, 1000000]
FRAGMENT_LENGTHS = [100000]
PATH = "/home/m4safari/projects/def-lila-ab/m4safari/ext2/data"


# PATH = "/content/drive/MyDrive/anew"


def experiment_task(args, env, fragment_length):
    print("\n Running the pipeline is started:")
    # Building the fragments
    fragment_file = f"{args['exp_type']}/fragments_{fragment_length}"
    print(f"\n Building fragment with length {fragment_length} is started.")
    run_fragment_builder(PATH, fragment_file, fragment_length, args['whole_genome'], env)
    print(f"\n Fragment with length {fragment_length} has been created.", flush=True)

    # Run the supervised classification under the first scenario (not challenging)
    result_folder = f"{PATH}/{args['exp_type']}/fragments_{fragment_length}"
    fasta_file = os.path.join(result_folder, env, f'Extremophiles_{env}.fas')
    print(f"\n deep clustering started.")
    run(NUM_CLUSTERS, fasta_file, args['max_k'], result_folder)
    print(f"\n deep clustering ended.", flush=True)



def run_pipeline(args):
    result = {}
    if args["exp_type"] == "exp1":
        print("clustering is running")
        for env in ENVS:
            for fragment_length in FRAGMENT_LENGTHS:
                result = experiment_task(args, env, fragment_length)


    # print(f"\n number of excuters: {os.cpu_count()}")
    # with ProcessPoolExecutor() as executor:
    #     futures = []
    #     for env in ENVS:
    #         for exp in range(num_exp):
    #             for fragment_length in FRAGMENT_LENGTHS:
    #                 # Submit tasks to the process pool
    #                 future = executor.submit(experiment_task, args, env, exp, fragment_length)
    #                 futures.append(future)
    #
    #     # Waiting for all futures to complete
    #     for future in as_completed(futures):
    #         future.result()  # You can handle exceptions here if needed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--exp_type', action='store', type=str)
    parser.add_argument('--max_k', action='store', type=int)
    parser.add_argument('--whole_genome', action='store_true')
    args = vars(parser.parse_args())

    run_pipeline(args)


if __name__ == '__main__':
    main()
