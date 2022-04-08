# mutable struct CallbackData
#     model::Gurobi.Optimizer
#     ptr::Ptr{Cvoid}
#     cb_where::Cint
# end

Base.cconvert(::Type{Ptr{Cvoid}}, x::CallbackData) = x
Base.unsafe_convert(::Type{Ptr{Cvoid}}, x::CallbackData) = x.ptr::Ptr{Cvoid}
Base.broadcastable(x::CallbackData) = Ref(x)

# mutable struct _CallbackUserData
#     model::Gurobi.Optimizer
#     callback::Function
# end
Base.cconvert(::Type{Ptr{Cvoid}}, x::_CallbackUserData) = x
function Base.unsafe_convert(::Type{Ptr{Cvoid}}, x::_CallbackUserData)
    return pointer_from_objref(x)::Ptr{Cvoid}
end

function _gurobi_callback_wrapper(
    p_model::Ptr{Cvoid},
    cb_data::Ptr{Cvoid},
    cb_where::Cint,
    p_user_data::Ptr{Cvoid}
)
    user_data = unsafe_pointer_to_objref(p_user_data)#::_CallbackUserData
    #Base.cconvert(_CallbackUserData, user_data)
    try
        user_data.callback(
            CallbackData(user_data.model, cb_data, cb_where),
            cb_where,
        )
    catch ex
        GRBterminate(p_model)
        if !(ex isa InterruptException)
            rethrow(ex)
        end
    end
    return Cint(0)
end

grb_callback = @cfunction(
    _gurobi_callback_wrapper,
    Cint,
    (Ptr{Cvoid}, Ptr{Cvoid}, Cint, Ptr{Cvoid})
)

function GRBtimedout(m::JuMP.Model)
    status = Ref{Cint}()
    GRBgetintattr(backend(m).optimizer.model, "Status", status)
    if status[] == 9
        return true
    else
        return false
    end
end
