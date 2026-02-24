module quasr
using PythonCall
using Downloads



function download_from_quasr(quasrID, tofile=quasrID)

    quasr_address = "https://quasr.flatironinstitute.org/simsopt_serials/"

    quasr_file = string(quasr_address, quasrID[1:4], "/", "serial", quasrID, ".json")

    to_file = string(pwd(), quasrID, ".json")

    quasrJSON = Downloads.download(
        quasr_file, to_file
    )
end

function load_from_quasr(quasrID)
    simsoptcore = pyimport("simsopt._core")

    quasr_file = string("serial", quasrID, ".json")

    simsoptcore.load(quasr_file)

    _, coils = simsoptcore.load(quasr_file)

    curve_fourier_amplitudes = [coils[1].curve.dofs_matrix[i] for i in 1:3]

    curve_fourier_names = coils[0].dof_names


end

end
