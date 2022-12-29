function [] = function_close_sutter( Sutter )

    fclose(Sutter.obj);

    clear('Sutter');
end

