from langchain_anthropic import ChatAnthropic
from langchain_openai import ChatOpenAI
from langchain.prompts import PromptTemplate
import pandas as pd
from os import listdir, path, makedirs
import subprocess
import json
import shutil
from typing import List, Tuple


LLMS = {
    "claude-opus": ChatAnthropic(model="claude-3-opus-20240229", temperature=0.0),
    "claude-sonnet": ChatAnthropic(model="claude-3-sonnet-20240229", temperature=0.0),
    "gpt-4o": ChatOpenAI(model="gpt-4o-2024-05-13", temperature=0.0),
    "gpt-turbo": ChatOpenAI(model="gpt-4-turbo-2024-04-09", temperature=0.0),
    "gpt-4": ChatOpenAI(model="gpt-4-0613", temperature=0.0),
}
DEFAULT_LLM = LLMS["claude-sonnet"]

bug_fix_prompt = (

)


def llm_bug_fix(script: str, llm=DEFAULT_LLM) -> str:
    """Use a LLM (currently defaults to Claude Sonnet for cost reasons) to provide a
    revised version of a possibly-buggy script."""
    prompt_template = PromptTemplate(
        input_variables=["script"],
        template=(
            "You are a helpful expert bioinformatician. "
            "I will provide you with a biological analysis script with zero or more bugs in it. "
            "You may assume that bugs are logical bioinformatics errors or programming bugs, but "
            "reading input files and serializing outputs are both correct. "
            "Fix any and all errors or bugs and respond with the corrected script. "
            "Your response should just be the corrected script, nothing else. Do not include "
            "any preamble, triple backticks code formatting, or anything else besides the corrected script.\n\n"
            "Example input:\n# compute sqrt of 8\nfrom math import sqrt\nsqtr(7)\n\n"
            "Example output:\n# compute sqrt of 8\nfrom math import sqrt\nsqrt(8)\n\n"
            "\n\nInput script:\n\n{script}"
        ),
    )
    sequence = prompt_template | llm
    return sequence.invoke(script).content


def _run_script(script_path, output_path):
    if script_path[-3:] == ".py":
        subprocess.run(["python", script_path, output_path])
    elif script_path[-2:].lower() == ".r":
        subprocess.run(
            ["R", "-f", script_path, "--args", output_path],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )

def test_modified_scripts(
        script_dir,
        correct_script="correct_script.py",
        overwrite_revised_script_dir=True,
        llm=DEFAULT_LLM,
    ):
    # create dir for LLM-fixed scripts
    revised_script_dir = path.join(script_dir, "llm_revised")
    if path.exists(revised_script_dir):
        if overwrite_revised_script_dir:
            # If overwrite_revised_scrpt_dir is True, remove the existing directory and its contents
            shutil.rmtree(revised_script_dir)
        else:
            # If overwrite_revised_scrpt_dir is False, throw an exception
            raise Exception(f"The directory '{revised_script_dir}' already exists. Set overwrite_revised_scrpt_dir to True to overwrite it.")
    makedirs(revised_script_dir, exist_ok=False)

    # run correct script to get output
    reference_output_path = path.join(script_dir, correct_script + ".outputjson")
    _run_script(path.join(script_dir, correct_script), reference_output_path)
    with open(reference_output_path, "r") as file_handle:
        reference_output = json.load(file_handle)

    # Iterate through each broken script.
    # For each, get the LLM-revised script and compare to reference results
    results = {}
    bad_scripts = [
        f
        for f in listdir(script_dir)
        if path.isfile(path.join(script_dir, f)) and not f.startswith(correct_script)
    ]
    for bad_script_name in bad_scripts:
        print(f"Working on script {bad_script_name}")
        try:
            # get the LLM-revised script and write to a file
            with open(path.join(script_dir, bad_script_name), "r") as file_handle:
                contents = file_handle.read()
            contents_llm_fixed = llm_bug_fix(contents, llm)
            revised_path = path.join(revised_script_dir, bad_script_name)
            with open(revised_path, "w+") as file_handle:
                file_handle.write(contents_llm_fixed)
            # run that revised script
            output_path = path.join(revised_script_dir, bad_script_name + ".outputjson")
            _run_script(revised_path, output_path)
        except Exception as e:
            print(e)
        # compare results to reference solution
        if path.exists(output_path):
            with open(output_path, "r") as file_handle:
                revised_output = json.load(file_handle)
            results[bad_script_name] = {
                "output_exists": True,
                "outputs_match": revised_output == reference_output,
            }
        else:
            results[bad_script_name] = {"output_exists": False, "outputs_match": False}
    return results

def compare_llms_one_script(script_dir, llms=LLMS, correct_script="correct_script.py"):
    results = {k: test_modified_scripts(script_dir, llm=v, correct_script=correct_script) for (k,v) in llms.items()}
    flattened_data = [
        {
            'llm': llm_model,
            'bug': test_case,
            'output_exists': result['output_exists'],
            'outputs_match': result['outputs_match']
        }
        for llm_model, test_cases in results.items()
        for test_case, result in test_cases.items()
    ]
    return pd.DataFrame(flattened_data)

def compare_llms_several_scripts(scripts: List[Tuple[str, str]], llms=LLMS):
    script_dfs = []
    for (script_dir, correct_script) in scripts:
        script_df = compare_llms_one_script(script_dir, llms=llms, correct_script=correct_script)
        script_df["dir"] = script_dir
        script_dfs.append(script_df)
    return pd.concat(script_dfs, ignore_index=True)
