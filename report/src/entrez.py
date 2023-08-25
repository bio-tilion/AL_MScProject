import pandas as pd
from Bio import Entrez
from datetime import datetime
import entrez_credentials
from time import sleep


# User credentials to pass to E-utilities
Entrez.email = entrez_credentials.email
Entrez.api_key = entrez_credentials.api_key


def my_elink(**kwargs):
    """
    Function to handle Elink function calls
    """
    handle = Entrez.elink(
        **kwargs
    )
    result = Entrez.read(handle)
    handle.close()
    return result


def get_citedby(pmid: str) -> int:
    """
    Function for returning the number of papers that cite the original input

    Argument
        pmid     Pubmed ID
    """
    param = {
        "dbfrom": "pubmed",
        "linkname": "pubmed_pubmed_citedin",
        "id": pmid,
    }

    result = my_elink(**param)
    out_citedby = len(result[0]["LinkSetDb"][0]["Link"] if result[0]["LinkSetDb"] else list())

    sleep(0.2)
    return out_citedby


def get_pub_date(pmid: str) -> str:
    """
    Function for returning the date an article has been published ("entrez" PubDate)

    Argument
        pmid     Pubmed ID
    """
    param = {
        "db": "pubmed",
        "id": pmid,
    }

    handle = Entrez.efetch(**param)
    result = Entrez.read(handle)
    handle.close()

    year = result["PubmedArticle"][0]["PubmedData"]["History"][0]["Year"]
    month = result["PubmedArticle"][0]["PubmedData"]["History"][0]["Month"]
    day = result["PubmedArticle"][0]["PubmedData"]["History"][0]["Day"]
    date = "-".join([year, month, day])

    sleep(0.2)
    return date


def get_linkDb_id(pmid: str, db: str, **kwargs) -> str | None:
    """
    Function for returning the ID in a different NCBI Db linked to the original paper (given a PMID)

    Arguments:
        pmid        Pubmed ID
        db          NCBI database

    Return a string with the new ID or None object if no IDs are found
    """
    param = {
        "dbfrom": "pubmed",
        "db": db,
        "id": pmid,
    }

    result = my_elink(**param, **kwargs)
    out_id = result[0]["LinkSetDb"][0]["Link"][0]["Id"] if result[0]["LinkSetDb"] else None

    sleep(0.2)
    return out_id


def get_pmid(doi: str) -> str:
    """
    Function for returning the PMID of a paper given its DOI

    Argument
        doi     DOI string
    """
    handle = Entrez.esearch(
        db="pubmed",
        term=doi,
        field="doi",
    )
    result = Entrez.read(handle)
    handle.close()

    sleep(0.2)
    return result["IdList"][0]


if __name__ == "__main__":
    table = pd.read_csv("report/src/oh2021.csv", sep=";")

    table["PMID"] = table["DOI"].apply(get_pmid)
    table["pub_date"] = table["PMID"].apply(get_pub_date)
    table["pub_date"] = pd.to_datetime(table["pub_date"], format="%Y-%m-%d")
    table["cited_by"] = table["PMID"].apply(get_citedby)
    table["avg_year_cit"] = table["cited_by"] / table["pub_date"].apply(lambda date: (datetime.today() - date).days / 365.25)
    table["geo_dataset_id"] = table["PMID"].apply(get_linkDb_id, db="gds")
    table["bioproject_id"] = table["PMID"].apply(get_linkDb_id, db="bioproject")

    table.to_csv("report/src/oh2021_metric.csv", index=False)
