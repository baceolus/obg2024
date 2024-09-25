from langchain_anthropic import ChatAnthropic
from langchain_openai import ChatOpenAI
from langchain.prompts import PromptTemplate
import pandas as pd
from os import listdir, path, makedirs, chdir
import subprocess
import json
import shutil
from typing import List, Tuple


LLMS = {
    "claude-opus": ChatAnthropic(model="claude-3-opus-20240229", temperature=0.0),
    "claude-sonnet": ChatAnthropic(model="claude-3-sonnet-20240229", temperature=0.0),
    "claude-35-sonnet": ChatAnthropic(model="claude-3-5-sonnet-20240620", temperature=0.0),
    "gpt-4o": ChatOpenAI(model="gpt-4o-2024-05-13", temperature=0.0),
    "gpt-turbo": ChatOpenAI(model="gpt-4-turbo-2024-04-09", temperature=0.0),
    "gpt-4": ChatOpenAI(model="gpt-4-0613", temperature=0.0),
}
DEFAULT_LLM = LLMS["claude-sonnet"]

prompt = """You are a helpful expert bioinformatician.
I will provide you with a task description. Respond with the source code for a
Python script that can be executed from a Linux command line to execute the task.
Your response should contain just the source code for the script, and nothing else.
In particular, do not include any preamble, triple backtick code formatting, or
anything else besides the source code. The only non-standard-library Python dependencies
your script may use are Numpy 1.26.4, Pandas 2.2.1, primer3 2.0.3, and Biopython 1.83.
Your script may call out to standard Linux (Ubuntu 22.04) shell commands as well as bowtie2
(version 2.5.4) and jellyfish (version 2.3.1).

Your response will be saved into bioscript.py which will be executed via `python bioscript.py`
using Python version 3.12.3.

Task:
{task}
"""


def llm_code_gen(task: str, llm=DEFAULT_LLM) -> str:
    """Use a LLM to construct a Python script for a given bioinformatics task."""
    prompt_template = PromptTemplate(
        input_variables=["task"],
        template=prompt
    )
    sequence = prompt_template | llm
    return sequence.invoke(task).content


def _run_script(script_path):
    subprocess.run(["python", script_path])

def challenge_llm(
        challenge_dir,
        overwrite_outputs=True,
        llm=DEFAULT_LLM,
        llm_name="LLM",
    ):
    """Challenge a LLM to produce a Python script to accomplish the challenge in challenge_dir,
    which must contain:
    * a task description in a file task.txt
    * a correct output in a file correct.output
    * input files; inputs and names should be described in the task text. All input filenames
        must start with "input_"
    
    This method executes several steps:
    1. create an output directory to hold the LLM's script (overwriting existing files
        if overwrite_outputs is True, otherwise raising an exception). Output dir name
        is determined by llm_name.
    2. use llm_code_gen to prompt the LLM for a script, and save it to the output
        directory as bioscript.py
    3. copy input files to the output directory
    4. change working directory to the output directory
    5. run the LLM-produced bioscript.py using _run_script
    6. change back to the working directory that was active when this method was called
    7. return a dictionary: {
        "output_exists": whether the LLM-produced bioscript.py ran and produced any output,
        "output_match": whether the produced output matches correct.output
    }"""
    # Step 1: Create output directory
    output_dir = path.join(challenge_dir, llm_name)
    if path.exists(output_dir):
        if overwrite_outputs:
            shutil.rmtree(output_dir)
        else:
            raise FileExistsError(f"Output directory {output_dir} already exists")
    makedirs(output_dir)

    # Step 2: Prompt LLM for script and save to output directory
    with open(path.join(challenge_dir, "task.txt"), "r") as f:
        task = f.read()
    script = llm_code_gen(task, llm)
    with open(path.join(output_dir, "bioscript.py"), "w") as f:
        f.write(script)

    # Step 3: Copy input files to output directory 
    for file in listdir(challenge_dir):
        if file.startswith("input_"):
            shutil.copy(path.join(challenge_dir, file), output_dir)

    # Steps 4-6: Run LLM script in output directory
    orig_dir = path.abspath(".")
    chdir(output_dir)
    _run_script("bioscript.py")
    chdir(orig_dir)

    # Step 7: Check if output exists and matches correct output
    output_file = path.join(output_dir, "output.txt")
    output_exists = path.exists(output_file)
    
    output_match = False
    if output_exists:
        with open(output_file, "r") as f:
            output = f.read().strip()
        with open(path.join(challenge_dir, "correct.output"), "r") as f:
            correct_output = f.read().strip()
        output_match = (output == correct_output)

    return {
        "output_exists": output_exists, 
        "output_match": output_match
    }

def compare_llms_one_challenge(challenge_dir, llms=LLMS, overwrite_outputs=True):
    results = {
        llm_name: challenge_llm(
            challenge_dir,
            overwrite_outputs=overwrite_outputs,
            llm=llm,
            llm_name=llm_name
        ) 
        for llm_name, llm in llms.items()
    }
    
    flattened_data = [
        {
            'challenge': challenge_dir,
            'llm': llm_name,
            'output_exists': result['output_exists'],
            'output_match': result['output_match']
        }
        for llm_name, result in results.items()
    ]
    
    return pd.DataFrame(flattened_data)

def compare_llms_several_challenges(challenge_dirs: List[str], llms=LLMS):
    challenge_dfs = []
    for challenge_dir in challenge_dirs:
        challenge_dfs.append(compare_llms_one_challenge(challenge_dir))
    return pd.concat(challenge_dfs, ignore_index=True)
