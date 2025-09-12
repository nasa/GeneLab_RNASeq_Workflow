#!/usr/bin/env python

import requests
import argparse
import re
import sys
import json

def get_glds_to_osd_mapping():
    """Get mapping from GLDS to OSD accessions using the search API."""
    search_url = "https://osdr.nasa.gov/osdr/data/search?ffield=Data+Source+Type&fvalue=cgene&size=5000"
    
    try:
        response = requests.get(search_url)
        response.raise_for_status()
        data = response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data from search API: {e}", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError:
        print("Error decoding JSON response from search API", file=sys.stderr)
        sys.exit(1)
    
    mapping = {}
    
    # Parse the search results to build GLDS -> OSD mapping
    studies = data.get('hits', {}).get('hits', [])
    for hit in studies:
        source = hit.get('_source', {})
        identifiers = source.get('Identifiers', '')
        osd_accession = source.get('Accession', '')
        
        # Extract GLDS accessions from identifiers
        glds_matches = re.findall(r'GLDS-\d+', identifiers)
        
        if glds_matches and osd_accession.startswith('OSD-'):
            for glds_accession in glds_matches:
                mapping[glds_accession] = osd_accession
    
    return mapping

def get_osd_and_glds(accession, api_url):
    osd_accession = None
    glds_accessions = []

    # Check if the accession is OSD or GLDS
    if accession.startswith('OSD-'):
        osd_accession = accession
        
        # For OSD input, try to find associated GLDS accessions from the mapping
        mapping = get_glds_to_osd_mapping()
        # Reverse lookup: find GLDS accessions that map to this OSD
        glds_accessions = [glds for glds, osd in mapping.items() if osd == accession]
        
        # If no GLDS found in mapping, try Biological Data API
        if not glds_accessions:
            try:
                response = requests.get(api_url)
                response.raise_for_status()
                data = response.json()
                
                # Search in wildcard endpoint results
                for osd_id, osd_data in data.items():
                    if osd_id == accession:
                        metadata = osd_data.get("metadata", {})
                        identifiers = metadata.get("identifiers", "")
                        glds_accessions = re.findall(r'GLDS-\d+', identifiers)
                        break
            except (requests.exceptions.RequestException, json.JSONDecodeError):
                pass  # Fall back to empty list
    
    elif accession.startswith('GLDS-'):
        glds_accessions = [accession]
        
        # Get GLDS -> OSD mapping from search API
        mapping = get_glds_to_osd_mapping()
        osd_accession = mapping.get(accession)
    else:
        print("Invalid accession format. Please use 'OSD-###' or 'GLDS-###'.", file=sys.stderr)
        sys.exit(1)

    if not osd_accession or not glds_accessions:
        print(f"No data found for {accession}", file=sys.stderr)
        sys.exit(1)

    return osd_accession, glds_accessions

def main():
    parser = argparse.ArgumentParser(description="Retrieve OSD and GLDS accessions.")
    parser.add_argument('--accession', required=True, help="Accession in the format 'OSD-###' or 'GLDS-###'")
    parser.add_argument('--api_url', required=True, help="OSDR API URL")
    args = parser.parse_args()

    osd_accession, glds_accessions = get_osd_and_glds(args.accession, args.api_url)

    # Output the results in a way that Nextflow can capture
    print(f"{osd_accession}")
    print(f"{','.join(glds_accessions)}")

if __name__ == "__main__":
    main()
