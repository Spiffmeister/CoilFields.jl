using CoilFields
using HTTP
using JSON

quasrID = 0009914

quasr_address = "https://quasr.flatironinstitute.org/simsopt_serials/0122/"
quasrJSON = HTTP.get(
    string(string(quasr_address, "serial"), string(quasrID, pad=7))
)
quasrstr = String(quasrJSON.body)
quasrdict = JSON.Parser.parse(quasrstr)
