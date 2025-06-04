{
  description = "Probabilistically Approximate DNF volume counter";
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";
  };
  outputs =
    {
      self,
      nixpkgs,
    }:
    let
      inherit (nixpkgs) lib;
      systems = lib.intersectLists lib.systems.flakeExposed lib.platforms.linux;
      forAllSystems = lib.genAttrs systems;
      nixpkgsFor = forAllSystems (system: nixpkgs.legacyPackages.${system});
      fs = lib.fileset;

      pepin-package =
        {
          stdenv,
          fetchFromGitHub,
          cmake,
          pkg-config,
          gmp,
          zlib,
        }:
        stdenv.mkDerivation {
          name = "pepin";
          src = fs.toSource {
            root = ./.;
            fileset = fs.unions [
              ./src
              ./CMakeLists.txt
              ./cmake
              ./scripts
              ./pepinConfig.cmake.in
            ];
          };

          nativeBuildInputs = [
            cmake
            pkg-config
          ];
          buildInputs = [
            gmp
            zlib
          ];
        };
    in
    {
      packages = forAllSystems (
        system:
        let
          pepin = nixpkgsFor.${system}.callPackage pepin-package {
          };
        in
        {
          inherit pepin;
          default = pepin;
        }
      );
    };
}
