# Containers

sc_pileup replaced the monolithic `Demuxafy.sif` with focused, pinned containers
(see `config.yaml` → `containers:`). Versions match what `Demuxafy.sif` shipped so
results are unchanged.

| logical name | image | source |
|---|---|---|
| `cellsnp`    | `quay.io/biocontainers/cellsnp-lite:1.2.2--hb0d9459_3` | public biocontainer |
| `samtools`   | `quay.io/biocontainers/samtools:1.15--h1170115_1`      | public biocontainer |
| `bcftools`   | `quay.io/biocontainers/bcftools:1.13--h3a49de5_0`       | public biocontainer |
| `souporcell` | `cumulusprod/souporcell:2.5`                            | Docker Hub (minimap2 + renamer/retag + samtools/bcftools/bedtools + vartrix) |
| `popscle`    | `ghcr.io/swarbricklab/sc_pileup-popscle:da70fc78`      | **built here** (`popscle/Dockerfile`) |

Everything except `popscle` is a public image referenced directly — no build step.

## Why popscle is built here

The only public popscle biocontainer is the older `0.1beta` release, whose
`dsc-pileup` writes the per-cell UMI table in a different format than the commit
`Demuxafy.sif` used (`statgen/popscle` master `@ da70fc78`). Building that exact
commit keeps `pileup.umi.gz` byte-identical to the atlas golden. bioconda has no
package for the commit, so it needs a Dockerfile (cf. imputation's
`containers/minimac4`).

## Building + pushing popscle

The image is published to GHCR under the `swarbricklab` org. Build and push with
a token that has `write:packages` (the from-source build needs a container
runtime — Docker/Podman locally, or absconda's remote builder):

```bash
# with Docker/Podman
docker build -t ghcr.io/swarbricklab/sc_pileup-popscle:da70fc78 containers/popscle
echo "$GHCR_TOKEN" | docker login ghcr.io -u <user> --password-stdin
docker push ghcr.io/swarbricklab/sc_pileup-popscle:da70fc78

# or via absconda's remote builder (offloads the build; see swarbricklab/absconda)
absconda build --template containers/popscle/Dockerfile \
  --repository ghcr.io/swarbricklab/sc_pileup-popscle --tag da70fc78 --push
```

Tag convention: the pinned upstream commit (`da70fc78`). A future "upgrade all
tool versions" pass would re-tag against newer popscle + biocontainers and
re-bless the golden.
