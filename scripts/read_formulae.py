#!/usr/bin/env python3
def get_formulas():
    with open("surge_input.txt", "r") as f:
        return [line.strip() for line in f if line.strip()]

if __name__ == "__main__":
    for formula in get_formulas():
        print(formula)
