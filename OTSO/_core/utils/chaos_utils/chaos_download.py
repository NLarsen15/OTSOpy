import requests
import json
import hashlib
from pathlib import Path


ZENODO_API = "https://zenodo.org/api"
TIMEOUT = 30


def md5_file(path):
    """Calculate the MD5 checksum of a local file."""
    md5 = hashlib.md5()

    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            md5.update(chunk)

    return md5.hexdigest()


def download_mat_only(doi, out_dir=None):
    """
    Download .mat files from a Zenodo record.

    Files are saved relative to this Python script unless
    an output directory is explicitly provided.

    Existing files are checked against Zenodo's checksum.
    Unchanged files are skipped.
    Changed files are downloaded again.

    Returns
    -------
    list[str]
        Paths to the downloaded/existing .mat files.
    """

    s = requests.Session()

    s.headers.update({
        "User-Agent": "CHAOS-mat-downloader-final"
    })

    # --------------------------------------------------
    # Resolve output directory
    # --------------------------------------------------
    # By default:
    #
    # project/
    #   script.py
    #   data/
    #
    # The data directory is therefore relative to the
    # location of this Python file, NOT the working directory.
    # --------------------------------------------------

    if out_dir is None:
        out_dir = Path(__file__).resolve().parent.parent.parent / "citation_generator"
    else:
        out_dir = Path(out_dir).resolve()

    out_dir.mkdir(parents=True, exist_ok=True)

    # --------------------------------------------------
    # Resolve DOI
    # --------------------------------------------------

    r = s.get(
        f"{ZENODO_API}/records",
        params={
            "q": f'doi:"{doi}" OR conceptdoi:"{doi}"',
            "size": 1
        },
        timeout=TIMEOUT
    )

    r.raise_for_status()

    hits = r.json().get("hits", {}).get("hits", [])

    if not hits:
        raise ValueError(f"DOI not found: {doi}")

    record_id = hits[0]["id"]

    # --------------------------------------------------
    # Fetch record metadata
    # --------------------------------------------------

    r = s.get(
        f"{ZENODO_API}/records/{record_id}",
        timeout=TIMEOUT
    )

    r.raise_for_status()

    record = r.json()

    metadata = record.get("metadata", {})

    # --------------------------------------------------
    # Print citation
    # --------------------------------------------------

    print()
    print("Chaos Model Metadata")
    print("--------------------")

    title = metadata.get("title", "Unknown title")

    record_doi = metadata.get("doi", doi)

    publication_date = metadata.get("publication_date", "")
    year = publication_date[:4] if publication_date else ""

    creator_names = []

    for creator in metadata.get("creators", []):

        if "person_or_org" in creator:
            person_or_org = creator["person_or_org"]

            if "name" in person_or_org:
                creator_names.append(person_or_org["name"])

        elif "name" in creator:
            creator_names.append(creator["name"])

    authors = ", ".join(creator_names)

    print(f"Title   : {title}")
    print(f"Authors : {authors}")
    print(f"DOI     : {record_doi}")
    print(f"Year    : {year}")

    # --------------------------------------------------
    # Save citation to JSON
    # --------------------------------------------------

    citation = (
        f"{authors} ({year}). "
        f"{title}. "
        f"Zenodo. https://doi.org/{record_doi}"
    )

    print()
    print("Citation:")
    print(citation)
    print()

    citation_path = out_dir / "citations.json"

    # Load existing citation database, or create it if it doesn't exist
    if citation_path.exists():
        with open(citation_path, "r", encoding="utf-8") as f:
            citation_data = json.load(f)
    else:
        citation_data = {"citations": {}}

    # Create a unique key for this citation
    citation_key = "CHAOS"

    # Append the new citation
    citation_data["citations"][citation_key] = {
        "citation": citation,
        "doi": record_doi
    }

    # Write the updated citation database back to the file
    with open(citation_path, "w", encoding="utf-8") as f:
        json.dump(
            citation_data,
            f,
            indent=2,
            ensure_ascii=False
        )

    print(f"Citation saved to: {citation_path}")



    # --------------------------------------------------
    # Find .mat files
    # --------------------------------------------------

    mat_files = []

    for file_info in record.get("files", []):

        name = file_info.get("key", "")

        if not name.lower().endswith(".mat"):
            continue

        url = file_info.get("links", {}).get("self")

        checksum = file_info.get("checksum", "")

        mat_files.append(
            (
                name,
                url,
                checksum
            )
        )

    if not mat_files:
        raise RuntimeError("No .mat files found in record.")

    print()
    print(f"Found {len(mat_files)} .mat file(s).")
    print()

    # --------------------------------------------------
    # Download / verify files
    # --------------------------------------------------

    downloaded = []

    datadir = Path(__file__).resolve().parent / "data"

    datadir.mkdir(parents=True, exist_ok=True)

    for name, url, zenodo_checksum in mat_files:

        path = datadir / name

        # --------------------------------------------------
        # Existing file
        # --------------------------------------------------

        if path.exists():

            # Zenodo normally provides checksums in the form:
            #
            # md5:xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            #
            if zenodo_checksum.startswith("md5:"):

                remote_hash = zenodo_checksum.split(":", 1)[1]

                local_hash = md5_file(path)

                if local_hash == remote_hash:

                    print(f"[same]    {name}")
                    print(f"          {path}")

                    downloaded.append(str(path))

                    continue

                else:

                    print(f"[changed] {name}")
                    print("          Existing file differs from Zenodo.")
                    print("          Downloading updated file...")

            else:

                print(f"[unknown] {name}")
                print("          Zenodo checksum unavailable.")
                print("          Downloading file...")

        # --------------------------------------------------
        # File does not exist
        # --------------------------------------------------

        else:

            print(f"[new]     {name}")
            print("          File does not exist locally.")
            print("          Downloading...")

        # --------------------------------------------------
        # Download
        # --------------------------------------------------

        try:

            with s.get(
                url,
                stream=True,
                timeout=TIMEOUT
            ) as r:

                r.raise_for_status()

                with open(path, "wb") as f:

                    for chunk in r.iter_content(
                        chunk_size=1024 * 1024
                    ):

                        if chunk:
                            f.write(chunk)

        except Exception:

            # Remove partially downloaded file if
            # something went wrong.
            if path.exists():
                path.unlink()

            raise

        # --------------------------------------------------
        # Verify downloaded file
        # --------------------------------------------------

        if zenodo_checksum.startswith("md5:"):

            remote_hash = zenodo_checksum.split(":", 1)[1]

            local_hash = md5_file(path)

            if local_hash != remote_hash:

                # Do not keep a corrupted/incomplete file.
                path.unlink()

                raise RuntimeError(
                    f"Checksum verification failed for {name}.\n"
                    f"Expected: {remote_hash}\n"
                    f"Received: {local_hash}"
                )

            print(f"[ok]      {name}")
            print("          Checksum verified.")

        else:

            print(f"[ok]      {name}")
            print("          Download complete.")

        downloaded.append(str(path))

        print()

    # --------------------------------------------------
    # Return downloaded/existing files
    # --------------------------------------------------

    return downloaded

