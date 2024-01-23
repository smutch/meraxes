{
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/612f97239e2cc474c13c9dafa0df378058c5ad8d";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        fftw-single = pkgs.fftw.override {
          enableMpi = false;
          precision = "single";
        };
        fftw-single-mpi = pkgs.fftw.override {
          enableMpi = true;
          precision = "single";
        };
      in
      with pkgs;
      {
        devShell = mkShell {
          buildInputs = [
            cmake
            openmpi
            gsl
            hdf5-mpi
            fftw-single  # <- required?
            fftw-single-mpi
            just
          ];
        };
      }
   );
}
