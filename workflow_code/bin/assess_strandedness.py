#! /usr/bin/env python

# a script to assess the dataset-wise strand selection from the sample-wise output from RSeQC infer_experiment
from pathlib import Path
from typing import Tuple
from statistics import median
import warnings
import argparse
import math

# DONE: Update this to be in sync with check configuration, as of 04/25/2022
STRANDEDNESS_ASSIGNMENT_THRESHOLD_MAX = 1.00 # values between above this will be assigned strandedness
STRANDEDNESS_ASSIGNMENT_THRESHOLD_MIN = 0.75
AMBIGUOUS_ASSIGNMENT_THRESHOLD_MAX = 0.75 # values above this, but below strandedness assignment will raise an warning
AMBIGUOUS_ASSIGNMENT_THRESHOLD_MIN = 0.60
UNSTRANDEDNESS_ASSIGMENT_THRESHOLD_MAX = 0.60 # values between this minimum and the ambiguous threshold are assigned unstranded
UNSTRANDEDNESS_ASSIGMENT_THRESHOLD_MIN = 0.40 # values between this minimum and the ambiguous threshold are assigned unstranded


def main(root_dir: str, strict=False):
    if strict:
        return assess_strict(root_dir)
    results = dict()
    sense_results = list()
    antisense_results = list()
    undetermined_results = list()
    skipped_files = []
    
    for file in Path(root_dir).glob("*"):
        with open(file, "r") as f:
            contents = f.read()
        
        # Skip files with Unknown Data type
        if "Unknown Data type" in contents:
            print(f"Warning: File {file.name} contains 'Unknown Data type' - excluding from strandedness assessment")
            skipped_files.append(file.name)
            continue
            
        result = _get_stranded_tuple(contents)
        undetermined_results.append(result[0])
        antisense_results.append(result[1])
        sense_results.append(result[2])
        results[file.name] = result
    
    # Check if we have any valid results
    if not sense_results:
        print("Warning: No valid strandedness data found. All samples had 'Unknown Data type' or invalid format.")
        # Default to unstranded if no valid data
        with open("result.txt", "w") as f:
            f.write("unstranded:0.5")
        return
    
    # determine average strandedness from valid samples only
    median_sense = median(sense_results)
    median_antisense = median(antisense_results)
    median_undetermined = median(undetermined_results)

    if median_sense > median_antisense:
        dominant, value = "sense", median_sense
    elif median_antisense >= median_sense:
        dominant, value = "antisense", median_antisense

    # assess strandedness
    if STRANDEDNESS_ASSIGNMENT_THRESHOLD_MAX > value > STRANDEDNESS_ASSIGNMENT_THRESHOLD_MIN:
        assignment = dominant
    elif AMBIGUOUS_ASSIGNMENT_THRESHOLD_MAX > value > AMBIGUOUS_ASSIGNMENT_THRESHOLD_MIN:
        warnings.warn(f"Strandedness assignment is not well defined for this dataset. median sense: {median_sense}, median antisense: {median_antisense}")
        assignment = dominant
    elif UNSTRANDEDNESS_ASSIGMENT_THRESHOLD_MAX > value > UNSTRANDEDNESS_ASSIGMENT_THRESHOLD_MIN:
        assignment = "unstranded"
    else:
        # Fallback in case the value doesn't fall into any of the defined ranges
        warnings.warn(f"Strandedness value ({value}) falls outside defined ranges. Defaulting to unstranded.")
        assignment = "unstranded"

    with open("result.txt", "w") as f:
        f.write(f"{assignment}:{value}")
    
    if skipped_files:
        print(f"Completed strandedness assessment. Skipped {len(skipped_files)} file(s) with 'Unknown Data type'.")


def _get_stranded_tuple(text: str) -> Tuple[float,float,float]:
    """ Parses stdout from infer_experiment """
    # Initialize with default values to avoid UnboundLocalError
    undetermined = 0.0
    sense = 0.0  
    antisense = 0.0
    
    for line in text.split("\n"):
        if line.startswith("Fraction of reads failed to determine:"):
            undetermined = float(line.split()[-1])
        elif line.startswith('Fraction of reads explained by "1++,1--,2+-,2-+":') or line.startswith('Fraction of reads explained by "++,--":') :
            antisense = float(line.split()[-1])
        elif line.startswith('Fraction of reads explained by "1+-,1-+,2++,2--":') or line.startswith('Fraction of reads explained by "+-,-+":'):
            sense = float(line.split()[-1])
    
    return (undetermined, sense, antisense)

def _strict_result(text, name):
    if "Unknown Data type" in text:
        raise ValueError(f"{name}: Unknown Data type; specify strandedness explicitly")
    prefixes = (
        ('Fraction of reads failed to determine:',),
        ('Fraction of reads explained by "1+-,1-+,2++,2--":',
         'Fraction of reads explained by "+-,-+":'),
        ('Fraction of reads explained by "1++,1--,2+-,2-+":',
         'Fraction of reads explained by "++,--":'),
    )
    values = []
    for alternatives in prefixes:
        lines = [line for line in text.splitlines() if line.startswith(alternatives)]
        if len(lines) != 1:
            raise ValueError(f"{name}: missing or duplicate inference fraction")
        try:
            value = float(lines[0].rsplit(':', 1)[1].strip())
        except ValueError:
            raise ValueError(f"{name}: malformed inference fraction") from None
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f"{name}: inference fractions must be finite and between 0 and 1")
        values.append(value)
    if not math.isclose(sum(values), 1.0, rel_tol=0, abs_tol=0.001):
        raise ValueError(f"{name}: inference fractions do not sum to 1")
    return tuple(values)


def _strict_assignment(sense, antisense):
    value = max(sense, antisense)
    if value >= STRANDEDNESS_ASSIGNMENT_THRESHOLD_MIN:
        return ('sense' if sense > antisense else 'antisense'), value
    if UNSTRANDEDNESS_ASSIGMENT_THRESHOLD_MIN <= value <= UNSTRANDEDNESS_ASSIGMENT_THRESHOLD_MAX:
        return 'unstranded', value
    raise ValueError(f"Ambiguous or insufficient strandedness inference ({value}); specify strandedness explicitly")


def assess_strict(root_dir):
    root = Path(root_dir)
    if not root.is_dir():
        raise ValueError(f"Inference directory does not exist: {root}")
    files = sorted(root.iterdir())
    if not files or any(not path.is_file() for path in files):
        raise ValueError("Expected one or more inference report files")
    results = [_strict_result(path.read_text(), path.name) for path in files]
    assignments = {_strict_assignment(row[2], row[1])[0] for row in results}
    if len(assignments) != 1:
        raise ValueError("Samples have inconsistent strandedness assignments; specify strandedness explicitly")
    assignment, value = _strict_assignment(median(row[2] for row in results), median(row[1] for row in results))
    Path('result.txt').write_text(f"{assignment}:{value}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('root_dir', nargs='?', default='infer_out')
    parser.add_argument('--strict', action='store_true', help='Reject unusable or ambiguous BAM inference reports')
    args = parser.parse_args()
    try:
        main(args.root_dir, strict=args.strict)
    except ValueError as error:
        if not args.strict:
            raise
        parser.exit(1, f"Strandedness assessment failed: {error}\n")
