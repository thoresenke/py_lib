"""
Playwright-only Scholar count helper.

This script opens a visible Chromium browser (headful) and navigates to the
Google Scholar query. If Google shows a CAPTCHA, the script will keep the
browser open and prompt you in the terminal to solve it; press Enter when
you've completed the CAPTCHA and the script will re-check and continue.

This approach intentionally uses an interactive, visible browser so you can
manually pass human verification when required.
"""

from typing import Optional, Dict
import re
import time
import random
import requests
import html as _html_lib
import os
from pathlib import Path
from playwright.sync_api import sync_playwright

def ScholarArticleCount2(
    search_string: str,
    year: int,
    timeout: int = 60,
    proxies: Optional[Dict[str, str]] = None,
    verbose: bool = True,
    user_data_dir: Optional[str] = None,
) -> Optional[int]:
    """Fetch Google Scholar result count using Playwright only.

    - Opens a visible browser so you can solve CAPTCHAs manually.
    - Waits for you to press Enter after you solve the CAPTCHA.
    """
    base_url = "https://scholar.google.com/scholar"
    params = {
        "q": search_string,
        "hl": "en",
        "as_sdt": "0,5",
        "as_ylo": str(year),
        "as_yhi": str(year),
    }

    headers = {
        "User-Agent": (
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/123.0.0.0 Safari/537.36"
        ),
    }

    number_re = re.compile(r"([0-9]+(?:[,\.\s\u00A0][0-9]{3})*)\s+results?", flags=re.IGNORECASE)
    block_markers = [
        "Our systems have detected unusual traffic",
        "To continue, please type the characters below",
        "detected unusual traffic from your network",
        "id=\"gs_captcha\"",
    ]

 
  
    # Default persistent profile directory (keeps cookies/localStorage between runs)
    if user_data_dir is None:
        repo_dir = Path(__file__).resolve().parent.parent
        user_data_dir = str(repo_dir / ".playwright_user_data")

    if verbose:
        print(f"Using Playwright user data dir: {user_data_dir}")

    os.makedirs(user_data_dir, exist_ok=True)

    with sync_playwright() as p:
        # Use a persistent context so CAPTCHA solutions and cookies are remembered
        context = p.chromium.launch_persistent_context(
            user_data_dir=user_data_dir,
            headless=False,
            user_agent=headers["User-Agent"],
            locale="en-US",
        )

        page = context.new_page()

        url = requests.Request("GET", base_url, params=params).prepare().url
        if verbose:
            print("Navigating to:", url)

        page.goto(url, timeout=timeout * 1000, wait_until="networkidle")

        # Loop until we find a count or the user confirms CAPTCHA cleared
        while True:
            time.sleep(0.5 + random.random() * 0.5)
            raw = page.content() or ""
            html_text = _html_lib.unescape(raw).replace('\u00A0', ' ').replace('\xa0', ' ')

            if verbose:
                snippet = html_text[:1000].replace('\n', ' ')
                print("HTML snippet (truncated):", snippet)

            # Detect blocking/CAPTCHA
            if any(marker in html_text for marker in block_markers):
                print("CAPTCHA or block detected in browser. Please solve it in the opened browser.")
                input("After solving the CAPTCHA, press Enter here to re-check the page...")
                # after user presses Enter, loop and re-evaluate
                continue

            # Try to extract number
            m = number_re.search(html_text)
            if m:
                count_str = re.sub(r"[,\.\s\u00A0]", "", m.group(1))
                try:
                    count = int(count_str)
                except ValueError:
                    count = None
                try:
                    context.close()
                except Exception:
                    pass
                return count

            # If no number and no CAPTCHA, ask the user whether to continue or quit
            if verbose:
                print("No result count found yet.")
            resp = input("Press Enter to retry page check, or type 'q' to quit: ")
            if resp.strip().lower() == 'q':
                try:
                    context.close()
                except Exception:
                    pass
                return None


if __name__ == "__main__":
    search_term = ' ""'
    years = list(range(2016, 2026))
    cnts =[]
    print("Open browser and solve CAPTCHA if present.\n")
    for year in years:
        cnt = ScholarArticleCount2(search_term, year, verbose=False)
        cnts.append(cnt) 
        print(f"Year {year}: {cnt}")
    
    print("\nSummary:")
    print("Year\tCount") 
    for year, cnt in zip(years, cnts):
        print(f"{year}\t{cnt}") 


