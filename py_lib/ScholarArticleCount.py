import requests
import html as _html_lib
import re
import time
import random
from typing import Optional, Dict


def ScholarArticleCount(
    search_string: str,
    year: int,
    pause: float = 2.0,
    max_retries: int = 3,
    proxies: Optional[Dict[str, str]] = None,
    use_playwright: bool = True,
    verbose: bool = True,
) -> Optional[int]:
    """Return an estimated result count from Google Scholar for a given year.

    Simple requests + retries and a Playwright fallback (if available).
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
            "Chrome/120.0.0.0 Safari/537.36"
        ),
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8",
        "Accept-Language": "en-US,en;q=0.9",
        "Referer": "https://www.google.com/",
    }

    session = requests.Session()

    number_re = re.compile(r"([0-9]+(?:[,\.\s\u00A0][0-9]{3})*)\s+results?", flags=re.IGNORECASE)
    block_markers = [
        "Our systems have detected unusual traffic",
        "To continue, please type the characters below",
        "detected unusual traffic from your network",
        "id=\"gs_captcha\"",
    ]

    # Try simple requests with retries/backoff
    for attempt in range(1, max_retries + 1):
        try:
            resp = session.get(base_url, params=params, headers=headers, timeout=15, proxies=proxies)
            if resp.status_code == 429:
                raise requests.HTTPError("429 Too Many Requests")
            resp.raise_for_status()
        except requests.RequestException as exc:
            if verbose:
                print(f"Request attempt {attempt} failed: {exc}")
            sleep_for = pause * (2 ** (attempt - 1)) + random.random()
            time.sleep(sleep_for)
            continue

        raw = resp.text or ""
        html_text = _html_lib.unescape(raw).replace('\u00A0', ' ').replace('\xa0', ' ')

        if verbose:
            print("HTTP status:", resp.status_code)
            print("HTML snippet:", html_text[:1000].replace("\n", " "))

        if any(marker in html_text for marker in block_markers):
            if verbose:
                print("Blocked or CAPTCHA detected in requests response")
            if attempt < max_retries:
                time.sleep(pause * (2 ** (attempt - 1)))
                continue
            break

        m = number_re.search(html_text)
        if not m:
            if verbose:
                print("Could not parse count from requests response")
            return None

        count_str = re.sub(r"[,\.\s\u00A0]", "", m.group(1))
        try:
            return int(count_str)
        except ValueError:
            return None

    # Playwright fallback
    if use_playwright:
        try:
            from playwright.sync_api import sync_playwright
        except Exception as exc:
            if verbose:
                print("Playwright not available:", exc)
            return None

        try:
            with sync_playwright() as p:
                browser = p.chromium.launch(headless=True)
                context = browser.new_context(user_agent=headers["User-Agent"], locale="en-US")
                page = context.new_page()
                page.goto(requests.Request("GET", base_url, params=params).prepare().url, timeout=30000, wait_until="networkidle")
                time.sleep(0.5 + random.random())
                raw = page.content() or ""
                browser.close()
        except Exception as exc:
            if verbose:
                print("Playwright navigation failed:", exc)
            return None

        html_text = _html_lib.unescape(raw).replace('\u00A0', ' ').replace('\xa0', ' ')
        if verbose:
            print("Playwright HTML snippet:", html_text[:1000].replace("\n", " "))

        if any(marker in html_text for marker in block_markers):
            if verbose:
                print("Blocked or CAPTCHA detected in Playwright response")
            return None

        m = number_re.search(html_text)
        if not m:
            if verbose:
                print("Could not parse count from Playwright response")
            return None

        count_str = re.sub(r"[,\.\s\u00A0]", "", m.group(1))
        try:
            return int(count_str)
        except ValueError:
            return None

    return None


if __name__ == "__main__":
    # Simple terminal output for years 2016-2025 for easy copy/paste
    search_term = '"fish"'
    years = range(2024, 2026)

    print("Year\tCount")
    for y in years:
        try:
            count = ScholarArticleCount(
                search_string=search_term,
                year=y,
                pause=1.6,
                max_retries=4,
                proxies=None,
                use_playwright=True,
                verbose=True,
            )
        # ...existing code...
        except Exception as exc: 
            count = None
            print(f"Error querying {y}: {exc}")

        print(f"{y}\t{'' if count is None else count}")
        time.sleep(0.5)