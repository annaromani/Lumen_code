#!/usr/bin/env python3
"""
Post-process FORD-generated documentation to inject custom sidebar navigation script
"""

import os
import re
from pathlib import Path

def inject_script(html_file):
    """Add custom sidebar script to HTML file if not already present"""
    with open(html_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Check if script is already injected
    if 'custom-sidebar.js' in content:
        return False
    
    # Check if this is a page file (documentation)
    if 'sidebar-toc' not in content:
        return False
    
    # Inject script before closing </head>
    script_tag = '    <script src="../js/custom-sidebar.js"></script>\n'
    if '</head>' in content:
        content = content.replace('</head>', f'{script_tag}</head>')
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(content)
        return True
    
    return False

def main():
    """Inject script into all generated page HTML files"""
    page_dir = Path('/Users/riccardo.reho/personalcodes/lumen/doc/ford_output/page')
    
    if not page_dir.exists():
        print(f"Error: Page directory not found: {page_dir}")
        return
    
    count = 0
    for html_file in page_dir.glob('*.html'):
        if inject_script(html_file):
            count += 1
            print(f"✓ Injected script into: {html_file.name}")
    
    print(f"\nTotal files updated: {count}")

if __name__ == '__main__':
    main()
