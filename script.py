import urllib.request, json
with urllib.request.urlopen("https://phasepro.elte.hu/download_full.json") as url:
    # Variable 'data' will contain the full database as a nested dictionary
    data = json.loads(url.read().decode())
    # Access the data of protein FUS
    print(data['P35637'])