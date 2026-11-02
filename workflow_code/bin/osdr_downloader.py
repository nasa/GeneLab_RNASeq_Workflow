#!/usr/bin/env python3

"""
OSDR File Downloader Script
Downloads files from NASA OSDR (Open Science Data Repository) based on specified parameters.

This script uses the OSDR biodata API to query and download files from GeneLab and ALSDA datasets.
It supports filtering by measurement type, technology type, and file extension.

Dependencies:
    - requests
    - argparse (built-in)
    - os, json, sys (built-in)

Version: 1.0
Author: Claude AI Assistant
Date: 2025-07-03
"""

import argparse
import requests
import json
import os
import sys
import urllib.parse
from pathlib import Path
import re
from typing import List, Dict, Optional, Tuple


class OSdRDownloader:
    """Class to handle OSDR file downloads."""
    
    def __init__(self):
        self.base_url = "https://visualization.osdr.nasa.gov/biodata/api/v2"
        self.session = requests.Session()
        
    def test_api_connectivity(self) -> bool:
        """Test if the API is accessible."""
        try:
            response = self.session.get(f"{self.base_url}/datasets/", timeout=10)
            response.raise_for_status()
            return True
        except requests.RequestException as e:
            print(f"Error: Cannot connect to OSDR API: {e}")
            return False
    
    def url_encode(self, text: str) -> str:
        """URL encode special characters."""
        return urllib.parse.quote(text.replace(' ', '%20'))
    
    def build_metadata_query_url(self, osd: str, measurement: Optional[str] = None, 
                                tech: Optional[str] = None, ext: Optional[str] = None, 
                                exclude_ext: Optional[str] = None, search: Optional[str] = None,
                                exclude_search: Optional[str] = None, genelab_only: bool = False) -> str:
        """Build the metadata query URL for the OSDR API."""
        
        # Start with basic required parameters for metadata endpoint
        query_parts = [
            f"id.accession={osd}",
            "file.file_name",
            "file.data_type", 
            "file.file_size",
            "file.remote_url",
            "file.category",  # Add category to identify GeneLab processed files
            "assay.protocol ref",  # Add protocol ref to check for GeneLab processing
            "format=json.records"  # Use proper JSON format for tables
        ]
        
        # Add GeneLab processed files filter if requested
        if genelab_only:
            query_parts.append("file.category=/GeneLab%20Processed%20.*%20Files/")
        
        # Add measurement type filter if specified
        if measurement:
            # Replace parentheses with wildcards for regex matching
            measurement_pattern = measurement.replace('(', '.*').replace(')', '.*')
            encoded_measurement = urllib.parse.quote(measurement_pattern)
            query_parts.append(f"investigation.study%20assays.study%20assay%20measurement%20type=/{encoded_measurement}/")
        
        # Add technology type filter if specified  
        if tech:
            # Replace parentheses with wildcards for regex matching
            tech_pattern = tech.replace('(', '.*').replace(')', '.*')
            encoded_tech = urllib.parse.quote(tech_pattern)
            query_parts.append(f"investigation.study%20assays.study%20assay%20technology%20type=/{encoded_tech}/")
        
        # Build query string manually to avoid encoding regex patterns
        encoded_parts = []
        regex_parts = []
        
        # First add all non-regex parts (these should be URL encoded)
        for part in query_parts:
            if not part.startswith('file.file_name=') and not part.startswith('file.file_name!='):
                encoded_parts.append(part)
        
        # Add file extension filter if specified (include only these extensions)
        if ext:
            regex_parts.append(f"file.file_name=/\\.{ext}$/")
        
        # Add exclude extension filter if specified (exclude these extensions)
        if exclude_ext:
            regex_parts.append(f"file.file_name!=/\\.{exclude_ext}$/")
        
        # Add search string filter if specified (filename contains this string)
        if search:
            regex_parts.append(f"file.file_name=/{search}/")
        
        # Add exclude search string filter if specified (filename does not contain this string)
        if exclude_search:
            regex_parts.append(f"file.file_name!={exclude_search}")
        
        # Combine all parts without encoding the regex patterns
        all_parts = encoded_parts + regex_parts
        query_string = '&'.join(all_parts)
        return f"{self.base_url}/query/metadata/?{query_string}"
    
    def build_file_download_url(self, filename: str, remote_url: str = "", osd: str = "") -> str:
        """Build URL to download a specific file using the remote_url from API response."""
        if remote_url:
            # Use the remote_url from API response and prepend the OSDR base URL
            if remote_url.startswith('/'):
                return f"https://osdr.nasa.gov{remote_url}"
            elif remote_url.startswith('http'):
                return remote_url
            else:
                return f"https://osdr.nasa.gov/{remote_url}"
        else:
            # Fallback to REST API if remote_url not available
            encoded_filename = urllib.parse.quote(filename)
            return f"{self.base_url}/query/data/?file.file_name={encoded_filename}"
    
    def is_genelab_processed(self, filename: str, data_type: str, file_category: str = "", protocol_ref: str = "") -> bool:
        """Check if a file is GeneLab processed data using protocol reference and file category."""
        
        # First check protocol reference (most reliable method for GeneLab processing)
        if protocol_ref:
            protocol_lower = protocol_ref.lower()
            if "genelab" in protocol_lower and "data processing protocol" in protocol_lower:
                return True
        
        # Second check file category
        if file_category:
            category_lower = file_category.lower()
            if "genelab processed" in category_lower and "files" in category_lower:
                return True
        
        # Fallback to filename/data type patterns if protocol and category not available
        # But be more conservative - only use obvious GeneLab patterns
        genelab_patterns = [
            r'^GLDS-\d+_.*_(unnormalized|normalized|differential).*counts',  # More specific count patterns
            r'^GLDS-\d+_.*_differential_expression',  # Differential expression files
            r'^GLDS-\d+_.*_(VST|RSEM|STAR)_.*counts',  # Specific analysis tool outputs
            r'^GLDS-\d+_.*_contrasts',  # Contrast files
            r'^GLDS-\d+_.*_sampletable',  # Sample tables (case insensitive)
        ]
        
        filename_lower = filename.lower()
        data_type_lower = data_type.lower() if data_type else ""
        
        # Check specific GeneLab filename patterns
        for pattern in genelab_patterns:
            if re.search(pattern, filename_lower):
                return True
        
        # Check data type for GeneLab-specific terms
        if data_type_lower:
            genelab_data_types = [
                'unnormalized counts',
                'normalized counts', 
                'differential expression',
                'sample table',
                'differential expression contrasts'
            ]
            if any(dt in data_type_lower for dt in genelab_data_types):
                return True
        
        return False
    
    def format_size(self, size_bytes: Optional[int]) -> str:
        """Format file size in human readable format."""
        if not size_bytes or size_bytes == 0:
            return "Unknown"
        
        for unit in ['B', 'KB', 'MB', 'GB']:
            if size_bytes < 1024:
                return f"{size_bytes:.1f}{unit}"
            size_bytes /= 1024
        return f"{size_bytes:.1f}TB"
    
    def get_measurement_tech_combinations(self, osd: str) -> List[Tuple[str, str]]:
        """Get all measurement/technology combinations for a dataset."""
        try:
            # Use metadata endpoint to discover available combinations
            url = f"{self.base_url}/query/metadata/?id.accession={osd}&investigation.study%20assays.study%20assay%20measurement%20type&investigation.study%20assays.study%20assay%20technology%20type&format=json.records"
            
            print(f"Discovering available data for {osd}...")
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            if isinstance(data, dict) and 'error' in data:
                print(f"API returned error: {data['error']}")
                return []
            
            if not data or len(data) == 0:
                print(f"No metadata found for {osd}")
                return []
            
            print(f"Found {len(data)} metadata records for {osd}")
            
            combinations = set()
            
            # Extract measurement/tech combinations from the records
            for record in data:
                measurement = record.get('investigation.study assays.study assay measurement type', '')
                tech = record.get('investigation.study assays.study assay technology type', '')
                
                if measurement and tech and measurement != 'null' and tech != 'null':
                    combinations.add((measurement, tech))
            
            # If no combinations found, try basic metadata query
            if not combinations:
                basic_url = f"{self.base_url}/query/metadata/?id.accession={osd}&file.file_name&file.data_type&format=json.records"
                response = self.session.get(basic_url, timeout=30)
                response.raise_for_status()
                basic_data = response.json()
                
                if basic_data:
                    print("Inferring measurement/tech combinations from available data...")
                    # Look for patterns in the data to infer measurement/tech types
                    for record in basic_data:
                        filename = record.get('file.file_name', '')
                        data_type = record.get('file.data_type', '')
                        
                        # Infer measurement/tech from filename patterns
                        if 'rna' in filename.lower() or 'rna' in data_type.lower():
                            if 'seq' in filename.lower() or 'seq' in data_type.lower():
                                combinations.add(("transcription profiling", "RNA sequencing"))
                            else:
                                combinations.add(("transcription profiling", "RNA-Seq"))
                        elif 'microarray' in filename.lower() or 'microarray' in data_type.lower():
                            combinations.add(("transcription profiling", "microarray"))
                        elif 'proteom' in filename.lower() or 'mass' in data_type.lower():
                            combinations.add(("protein expression profiling", "mass spectrometry"))
            
            if not combinations:
                # Default fallback - assume basic transcription profiling
                combinations.add(("transcription profiling", "RNA sequencing"))
            
            result = list(combinations)
            print(f"Discovered measurement/technology combinations: {result}")
            return result
            
        except Exception as e:
            print(f"Warning: Could not retrieve measurement/technology combinations: {e}")
            return [("transcription profiling", "RNA sequencing")]
    
    def query_files(self, osd: str, measurement: Optional[str] = None, 
                   tech: Optional[str] = None, ext: Optional[str] = None, 
                   exclude_ext: Optional[str] = None, search: Optional[str] = None,
                   exclude_search: Optional[str] = None) -> List[Dict]:
        """Query files from the OSDR API using metadata endpoint."""
        
        url = self.build_metadata_query_url(osd, measurement, tech, ext, exclude_ext, search, exclude_search)
        print(f"Querying: {url}")
        
        try:
            response = self.session.get(url, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            
            # Handle different response formats
            if isinstance(data, dict):
                if 'error' in data:
                    raise Exception(f"API Error: {data['error']}")
                # If it's a dict with a data field, extract the data
                if 'data' in data:
                    data = data['data']
                else:
                    raise Exception("Unexpected response format")
            
            if not isinstance(data, list):
                raise Exception(f"Expected list response, got {type(data)}")
            
            print(f"Found {len(data)} files")
            return data
            
        except requests.RequestException as e:
            raise Exception(f"Failed to query API: {e}")
    
    def filter_results_manually(self, data: List[Dict], measurement: str, tech: str, 
                               ext: Optional[str], exclude_ext: Optional[str],
                               search: Optional[str], exclude_search: Optional[str]) -> List[Dict]:
        """Manually filter results when additional filtering is needed."""
        if not data:
            return data
        
        filtered_data = []
        
        for record in data:
            filename = record.get('file.file_name', '')
            
            # Apply extension filter if specified (include only these)
            if ext and not filename.lower().endswith(f'.{ext.lower()}'):
                continue
            
            # Apply exclude extension filter if specified (exclude these)
            if exclude_ext and filename.lower().endswith(f'.{exclude_ext.lower()}'):
                continue
            
            # Apply search string filter if specified (filename contains this string)
            if search and search.lower() not in filename.lower():
                continue
            
            # Apply exclude search string filter if specified (filename does not contain this string)
            if exclude_search and exclude_search.lower() in filename.lower():
                continue
            
            # For measurement/tech filtering, we rely on the API query
            # since the records should already be filtered by the metadata endpoint
            filtered_data.append(record)
        
        filters_applied = []
        if ext:
            filters_applied.append(f"include .{ext}")
        if exclude_ext:
            filters_applied.append(f"exclude .{exclude_ext}")
        if search:
            filters_applied.append(f"include '{search}'")
        if exclude_search:
            filters_applied.append(f"exclude '{exclude_search}'")
        
        if filters_applied:
            filter_desc = " and ".join(filters_applied)
            print(f"Manually filtered {len(data)} files to {len(filtered_data)} files ({filter_desc})")
        
        return filtered_data
    
    def create_tsv_file(self, data: List[Dict], output_dir: str, osd: str, 
                       measurement: str = "unknown", tech: str = "unknown") -> str:
        """Create a TSV file with filename and download link information."""
        
        if not data:
            return ""
        
        # Create TSV filename
        subdir_name = f"{measurement}_{tech}".replace(" ", "_").replace("-", "_")
        tsv_filename = f"{osd}_{subdir_name}_file_list.tsv"
        tsv_path = os.path.join(output_dir, tsv_filename)
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        try:
            with open(tsv_path, 'w', encoding='utf-8') as f:
                # Write header
                f.write("Filename\tDownload_URL\tFile_Size\tData_Type\tGeneLab_Processed\n")
                
                # Process unique files (remove duplicates)
                seen_files = set()
                unique_files = []
                
                for record in data:
                    filename = record.get('file.file_name')
                    if filename and filename not in seen_files:
                        seen_files.add(filename)
                        unique_files.append(record)
                
                # Write file information
                for record in unique_files:
                    filename = record.get('file.file_name', '')
                    file_size = record.get('file.file_size', '')
                    data_type = record.get('file.data_type', '')
                    file_category = record.get('file.category', '')
                    protocol_ref = record.get('assay.protocol ref', '')
                    
                    if not filename:
                        continue
                    
                    # Build download URL using the remote_url from API response
                    remote_url = record.get('file.remote_url', '')
                    download_url = self.build_file_download_url(filename, remote_url, osd)
                    
                    # Check if GeneLab processed
                    is_genelab = self.is_genelab_processed(filename, data_type, file_category, protocol_ref)
                    genelab_status = "Yes" if is_genelab else "No"
                    
                    # Format file size
                    size_str = self.format_size(file_size) if file_size else "Unknown"
                    
                    # Write TSV row
                    f.write(f"{filename}\t{download_url}\t{size_str}\t{data_type}\t{genelab_status}\n")
            
            print(f"✓ Created TSV file: {tsv_path}")
            return tsv_path
            
        except Exception as e:
            print(f"✗ Failed to create TSV file: {e}")
            return ""
    
    def download_file(self, filename: str, filepath: str, file_record: Dict, osd: str = "") -> bool:
        """Download a single file using the remote_url from API response."""
        try:
            # Get remote_url from the file record
            remote_url = file_record.get('file.remote_url', '')
            download_url = self.build_file_download_url(filename, remote_url, osd)
            
            print(f"  Downloading {filename}...")
            print(f"  URL: {download_url}")
            
            response = self.session.get(download_url, stream=True, timeout=60)
            response.raise_for_status()
            
            # Create directory if it doesn't exist
            os.makedirs(os.path.dirname(filepath), exist_ok=True)
            
            with open(filepath, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            
            print(f"  ✓ Downloaded: {filename}")
            return True
            
        except Exception as e:
            print(f"  ✗ Failed to download {filename}: {e}")
            
            # Try alternative download using REST API if remote URL fails
            try:
                print(f"  Trying REST API fallback...")
                fallback_url = f"{self.base_url}/query/data/?file.file_name={urllib.parse.quote(filename)}"
                print(f"  Fallback URL: {fallback_url}")
                response = self.session.get(fallback_url, stream=True, timeout=60)
                response.raise_for_status()
                
                with open(filepath, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                
                print(f"  ✓ Downloaded via REST API fallback: {filename}")
                return True
                
            except Exception as e2:
                print(f"  ✗ REST API fallback also failed: {e2}")
                return False
    
    def process_files(self, data: List[Dict], output_dir: str, list_only: bool = False, 
                     measurement: str = "unknown", tech: str = "unknown", osd: str = "") -> Dict:
        """Process the files from API response."""
        
        if not data:
            print("No files found matching the criteria.")
            return {'total': 0, 'downloaded': 0, 'failed': 0, 'genelab': 0}
        
        # Create measurement_technology subdirectory
        subdir_name = f"{measurement}_{tech}".replace(" ", "_").replace("-", "_")
        measurement_tech_dir = os.path.join(output_dir, subdir_name)
        genelab_dir = os.path.join(measurement_tech_dir, "GeneLab_processed_data_files")
        
        # Remove duplicates based on filename
        seen_files = set()
        unique_files = []
        
        for record in data:
            filename = record.get('file.file_name')
            if filename and filename not in seen_files:
                seen_files.add(filename)
                unique_files.append(record)
        
        if len(data) != len(unique_files):
            print(f"Removed {len(data) - len(unique_files)} duplicate file entries")
        
        stats = {'total': 0, 'downloaded': 0, 'failed': 0, 'genelab': 0}
        
        print(f"\n{'='*60}")
        if list_only:
            print(f"Files available for {measurement} / {tech}:")
        else:
            print(f"Downloading files for {measurement} / {tech}:")
        print(f"{'='*60}")
        
        for record in unique_files:
            filename = record.get('file.file_name')
            file_size = record.get('file.file_size')
            data_type = record.get('file.data_type', '')
            file_category = record.get('file.category', '')
            protocol_ref = record.get('assay.protocol ref', '')
            remote_url = record.get('file.remote_url', '')
            
            if not filename:
                continue
            
            stats['total'] += 1
            
            # Check if GeneLab processed using protocol reference, file category, and patterns
            is_genelab = self.is_genelab_processed(filename, data_type, file_category, protocol_ref)
            if is_genelab:
                stats['genelab'] += 1
                file_marker = "[GeneLab]"
                target_dir = genelab_dir
            else:
                file_marker = "[Raw]    "
                target_dir = measurement_tech_dir
            
            size_str = self.format_size(file_size) if file_size else "Unknown"
            data_type_str = data_type if data_type else "Unknown"
            
            print(f"{file_marker} {filename} ({size_str}) - {data_type_str}")
            
            if not list_only:
                filepath = os.path.join(target_dir, filename)
                if self.download_file(filename, filepath, record, osd):
                    stats['downloaded'] += 1
                else:
                    stats['failed'] += 1
        
        # Create TSV file if in list mode
        if list_only and unique_files:
            self.create_tsv_file(unique_files, output_dir, osd, measurement, tech)
        
        return stats
    
    def run(self, osd: str, measurement: Optional[str] = None, tech: Optional[str] = None, 
            ext: Optional[str] = None, exclude_ext: Optional[str] = None, 
            search: Optional[str] = None, exclude_search: Optional[str] = None,
            output_dir: Optional[str] = None, list_only: bool = False):
        """Main execution method."""
        
        # Test API connectivity
        print("Testing API connectivity...")
        if not self.test_api_connectivity():
            sys.exit(1)
        print("✓ API connectivity test passed\n")
        
        # Set default output directory
        if not output_dir:
            output_dir = f"osdr_downloads_{osd}"
        
        total_stats = {'total': 0, 'downloaded': 0, 'failed': 0, 'genelab': 0}
        
        # If measurement and tech are specified, download for that specific combination
        if measurement and tech:
            try:
                data = self.query_files(osd, measurement, tech, ext, exclude_ext, search, exclude_search)
                stats = self.process_files(data, output_dir, list_only, measurement, tech, osd)
                
                # Update total stats
                for key in total_stats:
                    total_stats[key] += stats[key]
                    
            except Exception as e:
                print(f"Error processing {measurement}/{tech}: {e}")
        
        # If only measurement specified, get all tech types for that measurement
        elif measurement:
            combinations = self.get_measurement_tech_combinations(osd)
            matching_combinations = [(m, t) for m, t in combinations if m.lower() == measurement.lower()]
            
            if not matching_combinations:
                print(f"No technology types found for measurement '{measurement}' in {osd}")
                return
            
            for m, t in matching_combinations:
                try:
                    print(f"\nProcessing {m} / {t}...")
                    data = self.query_files(osd, m, t, ext, exclude_ext, search, exclude_search)
                    stats = self.process_files(data, output_dir, list_only, m, t, osd)
                    
                    # Update total stats
                    for key in total_stats:
                        total_stats[key] += stats[key]
                        
                except Exception as e:
                    print(f"Error processing {m}/{t}: {e}")
        
        # If neither specified, get all measurement/tech combinations
        else:
            combinations = self.get_measurement_tech_combinations(osd)
            
            if not combinations:
                print(f"No measurement/technology combinations found for {osd}")
                return
            
            for measurement, tech in combinations:
                try:
                    print(f"\nProcessing {measurement} / {tech}...")
                    data = self.query_files(osd, measurement, tech, ext, exclude_ext, search, exclude_search)
                    stats = self.process_files(data, output_dir, list_only, measurement, tech, osd)
                    
                    # Update total stats
                    for key in total_stats:
                        total_stats[key] += stats[key]
                        
                except Exception as e:
                    print(f"Error processing {measurement}/{tech}: {e}")
        
        # Print summary
        self.print_summary(osd, total_stats, output_dir, list_only, measurement, tech, ext, exclude_ext, search, exclude_search)
    
    def print_summary(self, osd: str, stats: Dict, output_dir: str, list_only: bool,
                     measurement: Optional[str], tech: Optional[str], ext: Optional[str], 
                     exclude_ext: Optional[str], search: Optional[str], exclude_search: Optional[str]):
        """Print summary statistics."""
        print(f"\n{'='*60}")
        print("Summary Report:")
        print(f"{'='*60}")
        print(f"Dataset: {osd}")
        print(f"Total files found: {stats['total']}")
        print(f"GeneLab processed files: {stats['genelab']}")
        print(f"Raw data files: {stats['total'] - stats['genelab']}")
        
        if not list_only:
            print(f"Files downloaded: {stats['downloaded']}")
            print(f"Failed downloads: {stats['failed']}")
            print(f"Output directory: {output_dir}")
            if stats['genelab'] > 0:
                print(f"GeneLab processed files saved to subdirectories: */GeneLab_processed_data_files")
        else:
            print(f"Output directory: {output_dir}")
            if stats['total'] > 0:
                print(f"TSV file(s) created with download URLs")
        
        print(f"\nApplied filters:")
        if measurement:
            print(f"  Measurement type: {measurement}")
        if tech:
            print(f"  Technology type: {tech}")
        if ext:
            print(f"  File extension (include): {ext}")
        if exclude_ext:
            print(f"  File extension (exclude): {exclude_ext}")
        if search:
            print(f"  Search string (include): {search}")
        if exclude_search:
            print(f"  Search string (exclude): {exclude_search}")
        if not measurement and not tech and not ext and not exclude_ext and not search and not exclude_search:
            print(f"  None (all files)")
        
        print(f"{'='*60}")


def main():
    """Main function to handle command line arguments."""
    parser = argparse.ArgumentParser(
        description="Download files from NASA OSDR (Open Science Data Repository)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # List all files for OSD-101
  python osdr_downloader.py --osd OSD-101 --list
  
  # Download transcription profiling RNA-Seq files
  python osdr_downloader.py --osd OSD-101 --measurement "transcription profiling" --tech "RNA-Seq"
  
  # Download only CSV files for a specific measurement type
  python osdr_downloader.py --osd OSD-101 --measurement "transcription profiling" --ext csv
  
  # Download all files except large tar.gz files
  python osdr_downloader.py --osd OSD-101 --exclude-ext tar.gz
  
  # Download files containing 'counts' in filename
  python osdr_downloader.py --osd OSD-101 --search counts
  
  # Download files excluding those with 'raw' in filename
  python osdr_downloader.py --osd OSD-101 --exclude-search raw
  
  # Download all files to custom directory
  python osdr_downloader.py --osd OSD-101 --out ./my_downloads
        """
    )
    
    parser.add_argument("--osd", required=True, 
                       help="OSD number (e.g., OSD-101)")
    parser.add_argument("--measurement", 
                       help="Measurement type (e.g., 'transcription profiling')")
    parser.add_argument("--tech", 
                       help="Technology type (e.g., 'RNA-Seq')")
    parser.add_argument("--ext", 
                       help="File extension to include (e.g., 'csv', 'fastq.gz')")
    parser.add_argument("--exclude-ext", 
                       help="File extension to exclude (e.g., 'tar.gz', 'zip')")
    parser.add_argument("--search", 
                       help="Search string to include in filename (e.g., 'counts', 'normalized')")
    parser.add_argument("--exclude-search", 
                       help="Search string to exclude from filename (e.g., 'raw', 'temp')")
    parser.add_argument("--out", 
                       help="Output directory (default: osdr_downloads_OSD-#)")
    parser.add_argument("--list", action="store_true", 
                       help="List files instead of downloading them")
    
    args = parser.parse_args()
    
    # Validate OSD format
    if not re.match(r'^OSD-\d+$', args.osd):
        print("Error: OSD must be in format OSD-XXX (e.g., OSD-101)")
        sys.exit(1)
    
    # Validate that ext and exclude_ext are not the same
    if args.ext and args.exclude_ext and args.ext.lower() == args.exclude_ext.lower():
        print("Error: --ext and --exclude-ext cannot be the same extension")
        sys.exit(1)
    
    # Validate that search and exclude_search are not the same
    if args.search and args.exclude_search and args.search.lower() == args.exclude_search.lower():
        print("Error: --search and --exclude-search cannot be the same string")
        sys.exit(1)
    
    # Create downloader and run
    downloader = OSdRDownloader()
    downloader.run(
        osd=args.osd,
        measurement=args.measurement,
        tech=args.tech,
        ext=args.ext,
        exclude_ext=args.exclude_ext,
        search=args.search,
        exclude_search=args.exclude_search,
        output_dir=args.out,
        list_only=args.list
    )


if __name__ == "__main__":
    main()
