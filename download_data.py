import os
import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor

csv_path = "/l/users/darya.taratynova/scbb_project/corrupted.csv"  
df = pd.read_csv(csv_path)

download_dir = "downloaded"
os.makedirs(download_dir, exist_ok=True)

def download_file(file_url):
    file_name = file_url.split("/")[-1]  
    file_path = os.path.join(download_dir, file_name)
    
    try:
        response = requests.get(file_url, stream=True)
        response.raise_for_status()
        with open(file_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"Downloaded: {file_path}")
    except requests.exceptions.RequestException as e:
        print(f"Failed to download {file_url}: {e}")

with ThreadPoolExecutor(max_workers=5) as executor: 
    executor.map(download_file, df["Download"])
