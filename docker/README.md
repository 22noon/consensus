# Running Allelescope with Docker

The Docker image contains Allelescope, pysam and samtools/htslib 1.19.2 on
Ubuntu 24.04. It is the supported way to run Allelescope on **Windows**, where
pysam and samtools are not available, and works the same way on macOS and
Linux.

The viewer itself runs in your normal web browser. The container builds the
browser instances and serves them; you open them at `http://localhost:8000/`.

## 1. Install Docker

- **Windows:** install [Docker Desktop](https://docs.docker.com/desktop/setup/install/windows-install/)
  and keep the default WSL 2 backend.
- **macOS:** install [Docker Desktop](https://docs.docker.com/desktop/setup/install/mac-install/).
- **Linux:** install [Docker Engine](https://docs.docker.com/engine/install/).

## 2. Build the image

From the repository root (the folder with `Dockerfile`):

```
docker build -t allelescope .
```

The first build downloads about 500 MB and takes a few minutes. Rebuild after
pulling a new version of the code.

## 3. Try the demo

```
docker run --rm -p 127.0.0.1:8000:8000 allelescope demo
```

Wait for `Demo ready`, then open <http://localhost:8000/demo/>. Press
`Ctrl+C` in the terminal to stop. See [`demo/README.md`](../demo/README.md) for
what to look at.

`-p 127.0.0.1:8000:8000` makes the viewer reachable only from your own
computer. To use another port, change the first number, for example
`-p 127.0.0.1:9000:8000`, and open `http://localhost:9000/demo/`.

## 4. Use your own data

The container sees one folder from your computer, mounted at `/data`. Put your
input files in it; the browser instances are written to `Browsers` inside it.

```
my-project/
├── inputs/
│   ├── sample.bam
│   ├── reference.fa
│   └── calls.vcf.gz
└── Browsers/            created by the build step
```

In the commands below, run them from inside `my-project`. All paths after the
image name are **container paths**, relative to `/data`, with forward slashes,
on every operating system.

### Build an instance

**Windows (PowerShell):**

```powershell
docker run --rm -v "${PWD}:/data" allelescope build -b inputs/sample.bam -r inputs/reference.fa -v inputs/calls.vcf.gz -o Browsers/sampleA -t 4
```

**Windows (Command Prompt):**

```bat
docker run --rm -v "%cd%:/data" allelescope build -b inputs/sample.bam -r inputs/reference.fa -v inputs/calls.vcf.gz -o Browsers/sampleA -t 4
```

**macOS / Linux:**

```bash
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/data" allelescope \
  build -b inputs/sample.bam -r inputs/reference.fa -v inputs/calls.vcf.gz \
        -o Browsers/sampleA -t 4
```

`build` takes the same options as `build_variant_browser` (see the
[main README](../allelescope/README.md)).

### Serve all instances

**Windows (PowerShell):** `docker run --rm -p 127.0.0.1:8000:8000 -v "${PWD}:/data" allelescope`

**Windows (Command Prompt):** `docker run --rm -p 127.0.0.1:8000:8000 -v "%cd%:/data" allelescope`

**macOS / Linux:**

```bash
docker run --rm --user "$(id -u):$(id -g)" -p 127.0.0.1:8000:8000 -v "$PWD:/data" allelescope
```

Then open `http://localhost:8000/sampleA/`. Every folder under `Browsers` is
served at its own address. Add `-d --name allelescope` to keep the server
running in the background and `docker stop allelescope` to stop it.

## Commands

| Command | What it does |
| --- | --- |
| *(none)* or `serve` | Serve every instance under `/data/Browsers` |
| `demo` | Build the bundled demo instance (first time only), then serve |
| `build ARGS` | Run `build_variant_browser ARGS` |
| anything else | Run that command, e.g. `bash`, `samtools`, `generate_viewer` |

Settings, passed with `-e NAME=value`:

| Variable | Default | Meaning |
| --- | --- | --- |
| `ALLELESCOPE_DATA_DIR` | `/data/Browsers` | Folder holding the instances |
| `ALLELESCOPE_PORT` | `8000` | Port inside the container (publish it with `-p`) |
| `ALLELESCOPE_THREADS` | `2` | samtools threads used by `demo` |

## Notes

- **Internet access in the browser:** the viewer page loads IGV.js and
  Bootstrap from `cdn.jsdelivr.net`, so the computer showing the page needs
  internet access, as without Docker.
- **File ownership on Linux:** without `--user`, files in `Browsers` are owned
  by user 1000. On Windows and macOS, Docker Desktop handles ownership and
  `--user` is not needed.
- **Speed on Windows:** files on a Windows drive (`C:\...`) are slower to read
  from a container than files inside WSL. For large BAMs, keep the project in
  your WSL home folder (`\\wsl$\Ubuntu\home\<you>\...`) and run the commands
  from a WSL terminal using the macOS / Linux form.
- **Instances are portable:** an instance folder built with Docker can be
  served by a normal installation and the other way round.
