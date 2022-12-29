classdef sutterController < handle

    properties
        sp
        reference
        name = 'sutterMP285'
        com = 'COM3';
        baud = 9600;
    end

    properties (Dependent)
        position
        velocity
    end

    properties (Constant)
        vscaleFactor = 10;
        stepMult = 25;
    end

    methods

        function obj = sutterController(com, baud)
            if nargin == 2
                obj.com = com;
                obj.baud = baud;
            elseif nargin == 1
                obj.baud = baud;  
            end
            
            fprintf('Initializing sutter... ')
            obj.sp = serialport(obj.com, obj.baud);
            setDTR(obj.sp, true); % set data terminal ready
            setRTS(obj.sp, true)
            configureTerminator(obj.sp, 'CR'); % set terminator to CR
            obj.sp.Timeout = 60; % long timeout allows large move commands to complete
            % sutter requires some voodoo when first turned on
            obj.velocity = 777;
            obj.velocity = 50;
            obj.updatePanel()
            fprintf('done.\r')
        end

        function getStatus(obj)
            % Sends the 'status' query to the Sutter and interprets the response.
            % Returns the current settings for step multiplier and velocity.
            fprintf('sutterMP285: get status info:\n')
            writeline(obj.sp, 's') % sends status command
            statusbytes = read(obj.sp, 32, "uint8");
            read(obj.sp, 1, 'int8'); % ignore CR

            % the value of STEP_MUL ("Multiplier yields msteps/nm") is at bytes 25 & 26
            stepMul=double(statusbytes(26))*256+double(statusbytes(25));

            % the value of "XSPEED"  and scale factor is at bytes 29 & 30
            if statusbytes(30) > 127
                vScaleFactor = 50;
            else
                vScaleFactor = 10;
            end

            currentVelocity=double(bitand(127,statusbytes(30)))*256+double(statusbytes(29)); % mm/s?

            fprintf(1,' step_mul (usteps/um):        %g\n',stepMul);
            fprintf(1,' velocity (usteps/sec):       %g\n',currentVelocity);
            fprintf(1,' vscale factor (usteps/step): %g\n',vScaleFactor);

        end

        function cmd = um2bytes(obj, xyz)
            xyz=xyz(:);
            xyz_bytes = typecast(int32(xyz .* obj.stepMult),'uint8')';
            cmd = [uint8(109) xyz_bytes uint8(13)];
        end

        function moveTo(obj, xyz)
            % move to raw location
            cmd = obj.um2bytes(xyz);
            obj.issueMove(cmd)
        end

        function moveAxis(obj, val, axis)
            % move axis generically
            movVec = [0 0 0];
            movVec(axis) = movVec(axis) + val;
            new_xyz = obj.reference + movVec;
            cmd = obj.um2bytes(new_xyz);
            obj.issueMove(cmd)
        end

        function moveZ(obj, z)
            % move in X rel to 
            new_xyz = obj.reference + [0 0 z];
            cmd = obj.um2bytes(new_xyz);
            obj.issueMove(cmd)
        end

        function moveToRef(obj)
            cmd = obj.um2bytes(obj.reference);
            obj.issueMove(cmd)
        end

        function issueMove(obj, cmd)
            write(obj.sp, cmd, 'uint8')
            cr = [];
            cr = read(obj.sp, 1, 'uint8');
            if isempty(cr)
                warning('Sutter did not finish moving before timeout!')
            end
        end

        function updatePanel(obj)
            % causes the Sutter to display the XYZ info on the front panel
            writeline(obj.sp, 'n'); % Sutter replies with a CR
            read(obj.sp, 1, 'int8'); % read and ignore the carriage return
        end

        function setOrigin(obj)
            writeline(obj.sp, 'o')
            read(obj.sp, 1, 'int8') % ignore CR
        end

        function setRef(obj)
            obj.reference = obj.position;
        end

        function set.velocity(obj, velocity)
            % Change velocity command 'V'xxCR where xx= unsigned short (16bit) int velocity
            % set by bits 14 to 0, and bit 15 indicates ustep resolution  0=10, 1=50 uSteps/step
            % V is ascii 86

            vv = typecast(uint16(velocity), 'uint8');
            cmd = [uint8(86) vv uint8(13)];
            write(obj.sp, cmd, 'uint8');

            % sutter replies with a CR, ignore it
            read(obj.sp, 1, 'int8');
        end

        function position = get.position(obj)
            % Read the current position of the Sutter MP-285
            % The Sutter sends binary data in units = 25usteps.um, i.e., *0.04 microns, & carriage return
            % (e.g. 'c' returns xxxxyyyyzzzz as 3 signed long ints, and ASCII 13 for CR)

            writeline(obj.sp, 'c');
            xyz = read(obj.sp, 3, 'int32');
            read(obj.sp, 1, 'int8'); % ignore CR
            position = xyz./obj.stepMult;
        end

    end

                
              

end