import sys, json, ssl
import os
from urllib import request
from urllib.error import HTTPError
from time import sleep

class InterProDownloader:
    def __init__(self):
        self.SILENT = False

    # this is an adapted version of the script provided by interpro at https://www.ebi.ac.uk/interpro/result/download/#/protein/UniProt/entry/InterPro/IPR005527/|fasta
    def download_entry_proteins(self, interpro_entry_id, path):
        '''

        :param interpro_entry_id:
        :param path:
        :return:
        '''

        # check if file exists already
        if os.path.isfile(path):
            print('File already downloaded. EXIT. ')
            return

        base_url = 'https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/' + str(interpro_entry_id) + '/?page_size=200&extra_fields=sequence'
        if not self.SILENT:
            print('Downloading Data. This might take a while...')
        self._output_list(base_url, path)
        if not self.SILENT:
            print('Data successfully downloaded to', path)


    def _output_list(self, base_url, path):
        '''

        :param base_url:
        :param path:
        :return:
        '''

        HEADER_SEPARATOR = "|"
        LINE_LENGTH = 80

        #disable SSL verification to avoid config issues
        context = ssl._create_unverified_context()

        next = base_url
        last_page = False

        attempts = 0
        while next:
            try:
                req = request.Request(next, headers={"Accept": "application/json"})
                res = request.urlopen(req, context=context)
                # If the API times out due a long running query
                if res.status == 408:
                    # wait just over a minute
                    sleep(61)
                    # then continue this loop with the same URL
                    continue
                elif res.status == 204:
                    #no data so leave loop
                    break
                payload = json.loads(res.read().decode())
                next = payload["next"]
                attempts = 0
                if not next:
                    last_page = True
            except HTTPError as e:
                if e.code == 408:
                    sleep(61)
                    continue
                else:
                    # If there is a different HTTP error, it wil re-try 3 times before failing
                    if attempts < 3:
                        attempts += 1
                        sleep(61)
                        continue
                    else:
                        sys.stderr.write("LAST URL: " + next)
                        raise e

            for i, item in enumerate(payload["results"]):

                if ("entries" in item):
                    for entry in item["entries"]:
                        for locations in entry["entry_protein_locations"]:
                            for fragment in locations["fragments"]:
                                start = fragment["start"]
                                end = fragment["end"]


                                open(path, 'a').write(">" + item["metadata"]["accession"] + HEADER_SEPARATOR
                                                 + entry["accession"] + HEADER_SEPARATOR
                                                 + str(start) + "..." + str(end) + HEADER_SEPARATOR
                                                 + item["metadata"]["name"] + "\n")
                                seq = item["extra_fields"]["sequence"]
                                fastaSeqFragments = [seq[0+i:LINE_LENGTH+i] for i in range(0, len(seq), LINE_LENGTH)]
                                for fastaSeqFragment in fastaSeqFragments:
                                    open(path, 'a').write(fastaSeqFragment + "\n")
                else:
                    open(path, 'a').write(">" + item["metadata"]["accession"] + HEADER_SEPARATOR + item["metadata"]["name"] + "\n")
                    seq = item["extra_fields"]["sequence"]
                    fastaSeqFragments = [seq[0+i:LINE_LENGTH+i] for i in range(0, len(seq), LINE_LENGTH)]
                    for fastaSeqFragment in fastaSeqFragments:
                        open(path, 'a').write(fastaSeqFragment + "\n")

                # Don't overload the server, give it time before asking for more
            if next:
                sleep(1)
