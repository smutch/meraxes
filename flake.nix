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
        compileInputs = with pkgs; [
            openmpi
            gsl
            hdf5-mpi
            fftw-single-mpi
        ];
      in
      with pkgs;
      {
        devShell = mkShell {
          buildInputs = [
            cmake
            fftw-single  # <- required?
            just
          ] ++ compileInputs;

          clangdConf = pkgs.writeText ".clangd" (
            ''
            CompileFlags:
              Add: 
            '' + (
              builtins.concatStringsSep "\n" (
                map ( x: "    - \"--include-directory=" + ( lib.makeSearchPathOutput "dev" "include" [x] ) + "\"" ) compileInputs
              )
            )
          );

          shellHook = ''
            cp $clangdConf .clangd
            '';
        };

      }
   );
}
