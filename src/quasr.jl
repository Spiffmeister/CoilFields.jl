
function setup_for_pythoncall()
    ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
    ENV["JULIA_PYTHONCALL_EXE"] = string(pwd(), "/.venv/bin/python")
end

function download_from_quasr end
function load_from_quasr end
