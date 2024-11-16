{
  inputs = {
    nixpkgs = {
      url = "github:nixos/nixpkgs/nixos-unstable";
    };
    flake-utils = {
      url = "github:numtide/flake-utils";
    };
  };

  outputs = { nixpkgs, flake-utils, ... }: flake-utils.lib.eachDefaultSystem (system:
    let
      pkgs = import nixpkgs {
        inherit system;
      };
    in rec {
      devShell = pkgs.mkShell {
        buildInputs = with pkgs; [
          (python312.withPackages(ps: with ps; [
            ipython
            ipykernel
            bash-kernel
            snakemake
            bwa
            samtools
            pysam
            graphviz
            networkx
            matplotlib
            numpy
            rdkit
            jinja2
	    pip
          ]))
          jupyter
        ];

      };
    }
  );
}
